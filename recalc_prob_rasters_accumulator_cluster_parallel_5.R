library(tidyverse)
library(terra)
library(tidyterra)
library(optparse)
library(furrr)
library(future)
#library(pryr) #Install the package if you want to use the mem_used() function

#Set up input arguments with optparse
option_list = list(
  make_option(c("-w", "--working_directory"), type="character", default=NULL,
              help="working directory (mandatory)", metavar="character"),
  make_option(c("-i", "--foa_lcp_path"), type="character", default=NULL,
              help="foa lcp path (mandatory)", metavar="character"),
  make_option(c("-r", "--foa_run"), type="character", default=NULL,
              help="foa run label (mandatory)", metavar="character"),
  make_option(c("-s", "--scenario"), type="character", default=NULL,
              help="project scenario (mandatory)", metavar="character"),
  make_option(c("-t", "--run_timepoint"), type="character", default=NULL,
              help="timepoint for the scenario (mandatory)", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#Set working directory
setwd(opt$working_directory)
wd <- getwd()

#Set up temporary directory
temp_dir <- file.path(wd, "temp_dir")
dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(TMPDIR = temp_dir)

#Force Terra to do chunked processing
#terra::terraOptions(memfrac = 0.5, tempdir = temp_dir)  # Use only 50% of available memory

# Specify log file path
log_file <- paste0("./recalc_probability_rasters_", opt$foa_run, "_", opt$scenario, ".log")

# Helper function to log messages to the log file
log_message <- function(message) {
  cat(message, "\n", file = log_file, append = TRUE)
}

# Helper function to clean up temp directory 
cleanup_tempdir <- function(temp_dir) {
  if (dir.exists(temp_dir)) {
    log_message(paste0("Cleaning up temporary directory: ", temp_dir))
    tryCatch({
      unlink(temp_dir, recursive = TRUE, force = TRUE)
      log_message("Temporary directory successfully deleted.")
    }, error = function(e) {
      log_message(paste0("Error while deleting temp directory: ", e$message))
    })
  } else {
    log_message("Temp directory does not exist or was already cleaned up.")
  }
}


#Define a function to calc prob rasters using an accumulator method
#Note: accum_bp is a single-band accumulator raster. fl_accumulators is a list of six single-band accumulator rasters
calc_prob_w_accumulator <- function(season_fire_path, categories) {
  library(terra)
  season_id <- stringr::str_extract(season_fire_path, "(?<=Season)\\d+(?=_)")
  log_message(paste0("Now processing season ", season_id, "..."))
  seasonfire_FLs <- terra::rast(season_fire_path, lyr = 3)
  burned_mask <- !is.na(seasonfire_FLs)
  accum_bp <- burned_mask
  #seasonfire_FLs_int <- floor(seasonfire_FLs)
  accum_bp_path <- file.path(temp_dir, paste0("season", season_id, "_accum_bp.tif"))
  #log_message(print(accum_bp_path))
  terra::writeRaster(accum_bp, accum_bp_path, overwrite=TRUE, datatype="INT1U")
  
  fl_paths <- map2(categories, names(categories), function(bounds, name) {
    mask <- terra::ifel(!is.na(seasonfire_FLs) &
                          seasonfire_FLs >= bounds[1] &
                          seasonfire_FLs < bounds[2], TRUE, FALSE)
    log_message(paste0("Category ", name, ": ", global(mask, "sum", na.rm=TRUE)[[1]], " matching pixels"))
    acc <- terra::ifel(mask, 1, NA)
    path <- file.path(temp_dir, paste0("season", season_id, "_accum_fl_", name, ".tif"))
    terra::writeRaster(acc, path, overwrite=TRUE)
    return(path)
  }) |> set_names(names(categories))
  
  return(c(list(accum_bp = accum_bp_path), fl_paths))
}

#List the SeasonFire raster paths
season_fire_files <- list.files(path = paste0(wd, "/SeasonFires_merged_tifs_", opt$scenario),
                                pattern = ".tif$", full.names=TRUE)
r <- terra::rast(season_fire_files[1])
num_seasons <- length(season_fire_files)

#Define flame length categories
categories <- list(
  "0_2"    = c(0, 0.6096),
  "2_4"    = c(0.6096, 1.2192),
  "4_6"    = c(1.2192, 1.8288),
  "6_8"    = c(1.8288, 2.4384),
  "8_12"   = c(2.4384, 3.6576),
  "12plus" = c(3.6576, Inf)
)

# Run each seasonfire file in parallel
log_message("Calculating accumulator bp and flp rasters in parallel...")
n_workers <- 60
plan(cluster, workers = n_workers)
log_message(paste0("Launching with ", n_workers, " workers using PSOCK cluster..."))
start_time <- Sys.time()
log_message(paste0("Script started at: ", start_time))
#Set up global future options
furrr_options <- furrr_options(globals=c("wd", "opt", "categories", "season_fire_files", "num_seasons", 
                                         "calc_prob_w_accumulator", "log_message", "log_file", "temp_dir"),
                               packages=c("terra","tidyverse","tidyterra","stringr"), seed=TRUE)
results_list <- future_map(season_fire_files, ~calc_prob_w_accumulator(.x, categories),
                           .options=furrr_options,
                           .progress = FALSE)
#How much memory was used?
#log_message(paste0("Memory used (Mb): ", mem_used()/1024^2))
plan(sequential)
# Convert paths to SpatRaster objects
log_message("Reading temporary accumulator rasters from disk and combining...")
timepoint1 <- Sys.time()
log_message(paste0("System time before combining bp accumulator rasters: ", timepoint1))
# Calculate and log duration
duration1 <- difftime(timepoint1, start_time, units = "mins")
log_message(paste0("Duration: ", round(duration1, 2), " minutes"))

# Create accumulator rasters on disk, initialize with zeros
#template <- terra::rast(opt$foa_lcp_path, lyrs=1)
#terra::crs(template) <- terra::crs(template)
#terra::values(template) <- 0

plan(multisession, workers = 60)  # Adjust as needed

result_chunks <- split(results_list, ceiling(seq_along(results_list) / 2000))

partial_sums <- future_map(result_chunks, function(chunk) {
  acc_bp <- terra::setValues(terra::rast(terra::rast(opt$foa_lcp_path, lyrs=1)), 0)
  acc_flp <- map(names(categories), ~ {
    terra::setValues(terra::rast(opt$foa_lcp_path, lyrs=1), 0)
  }) |> set_names(names(categories))
  
  for (res in chunk) {
    # --- Burn probability raster ---
    season_bp <- terra::rast(res$accum_bp)
    acc_bp <- sum(acc_bp, season_bp, na.rm = TRUE)
    # --- Flame length category rasters ---
    for (cat in names(categories)) {
      season_fl <- terra::rast(res[[cat]])
      acc_flp[[cat]] <- sum(acc_flp[[cat]], season_fl, na.rm = TRUE)
    }
  }
  # Return both types
  list(bp = acc_bp, flp = acc_flp)
}, .progress = TRUE)


timepoint2 <- Sys.time()
log_message(paste0("System time after accumulator rasters and before writing: ", timepoint2))
# Calculate and log duration
duration2 <- difftime(timepoint2, start_time, units = "mins")
log_message(paste0("Duration: ", round(duration2, 2), " minutes"))

# Combine all the partial sums
log_message("Combining final_acc_bp...")
final_acc_bp <- reduce(map(partial_sums, "bp"), `+`)
log_message("Combining final_acc_flp...")
final_acc_flp <- reduce(map(partial_sums, "flp"), function(a, b) {
  map2(a, b, `+`)
})

log_message("Assigning final_acc_bp to accum_bp...")
accum_bp <- final_acc_bp
names(accum_bp) <- "recalc_bp"

# Calculate and write burn probability
log_message("Calculating burn probability raster...")
burn_prob <- accum_bp / num_seasons
log_message("Writing burn probability raster to disk...")
terra::writeRaster(burn_prob, filename=paste0("./recalc_bp_", opt$foa_run, "_", opt$scenario, "_", opt$run_timepoint, ".tif"), overwrite=TRUE, datatype="FLT4S")
log_message("Burn probability raster written.")

# Final FLP outputs
for (cat in names(categories)) {
  accum_fl <- final_acc_flp[[cat]]
  cflp <- accum_fl / accum_bp
  flp <- cflp * burn_prob
  log_message(paste0("Writing conditional flame length probability raster for flame lengths ", cat))
  terra::writeRaster(cflp, filename=paste0("./recalc_cflp_", cat, "_", opt$foa_run, "_", opt$scenario, "_", opt$run_timepoint, ".tif"), overwrite=TRUE, datatype="FLT4S")
  log_message(paste0("Writing unconditional flame length probability raster for flame lengths ", cat))
  terra::writeRaster(flp, filename=paste0("./recalc_flp_", cat, "_", opt$foa_run, "_", opt$scenario, "_", opt$run_timepoint, ".tif"), overwrite=TRUE, datatype="FLT4S")
}

end_time <- Sys.time()
log_message(paste0("Script ended at: ", end_time))

# Calculate and log duration
duration3 <- difftime(end_time, start_time, units = "mins")
log_message(paste0("Total duration: ", round(duration3, 2), " minutes"))

# Clean up parallel backend
plan(sequential)
# Clean up temporary directory
cleanup_tempdir(temp_dir)
