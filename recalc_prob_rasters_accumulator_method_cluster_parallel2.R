library(tidyverse)
library(terra)
library(tidyterra)
library(optparse)
library(furrr)
library(future)

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

# Specify log file path
log_file <- paste0("./recalc_probability_rasters_", opt$foa_run, "_", opt$scenario, ".log")

# Helper function to log messages to the log file
log_message <- function(message) {
  cat(message, "\n", file = log_file, append = TRUE)
}

#Define a function to calc prob rasters using an accumulator method
#Note: accum_bp is a single-band accumulator raster. fl_accumulators is a list of six single-band accumulator rasters
calc_prob_w_accumulator <- function(season_fire_path, categories, foa_lcp_path) {
  foa_lcp <- terra::rast(opt$foa_lcp_path, lyrs=1)
  foa_lcp <- terra::unwrap(foa_lcp)
  season_id <- stringr::str_extract(season_fire_path, "(?<=Season)\\d+(?=_)")
  log_message(paste0("Now processing season ", season_id, "..."))
  seasonfire_FLs <- terra::rast(season_fire_path, lyr = 3)
  terra::crs(seasonfire_FLs) <- terra::crs(foa_lcp)
  seasonfire_FLs <- terra::extend(seasonfire_FLs, terra::ext(foa_lcp), snap = "near")
  #log_message(paste0("FL raster summary: ", summary(seasonfire_FLs)))
  #log_message(paste0("FL min/max: ", terra::minmax(seasonfire_FLs)))
  burned_mask <- !is.na(seasonfire_FLs)
  #log_message(paste0("burned_mask: ", sum(burned_mask[], na.rm=TRUE), " burned pixels"))
  accum_bp <- burned_mask
  seasonfire_FLs_int <- floor(seasonfire_FLs)

  temp_dir <- file.path(tempdir(), paste0("accum_", tools::file_path_sans_ext(basename(season_fire_path))))
  dir.create(temp_dir, showWarnings = FALSE)
  
  accum_bp_path <- file.path(temp_dir, paste0("season", season_id, "_accum_bp.tif"))
  log_message(print(accum_bp_path))
  terra::writeRaster(accum_bp, accum_bp_path, overwrite=TRUE)

  fl_paths <- map2(categories, names(categories), function(bounds, name) {
    mask <- seasonfire_FLs_int >= bounds[1] & seasonfire_FLs_int < bounds[2]
    acc <- terra::ifel(mask, 1, NA)
    path <- file.path(temp_dir, paste0("season", season_id, "_accum_fl_", name, ".tif"))
    log_message(print(path))
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
  "0_2"    = c(0, 2),
  "2_4"    = c(2, 4),
  "4_6"    = c(4, 6),
  "6_8"    = c(6, 8),
  "8_12"   = c(8, 12),
  "12plus" = c(12, Inf)
)

# Run each seasonfire file in parallel
log_message("Calculating accumulator bp and flp rasters in parallel...")
plan(multisession)
#Set up global future options
furrr_options <- furrr_options(globals=c("wd", "opt", "categories", "season_fire_files", "num_seasons", "calc_prob_w_accumulator", "log_message", "log_file"), seed=TRUE)
results_list <- future_map(season_fire_files, ~calc_prob_w_accumulator(.x, categories, foa_lcp_path),
                           .options=furrr_options,
                           .progress = TRUE)

# Convert paths to SpatRaster objects
log_message("Reading temporary accumulator rasters from disk and combining...")

# Convert paths into SpatRasters
results_rasters <- map(results_list, function(res) {
  map(res, rast)  # each res is a list of paths like accum_bp and fl categories
})

combined_results <- Reduce(function(x, y) {
  Map(function(a, b) {
    terra::ifel(!is.na(a), a + b, b)  # safely add only non-NA cells
  }, x, y)
}, results_rasters)

# Unpack combined results
accum_bp <- combined_results$accum_bp
accum_flps <- combined_results[names(categories)]
names(accum_bp) <- "recalc_bp"

#Calculate the burn probability raster
burn_prob <- accum_bp/num_seasons
#Write the burn probability raster
log_message("Writing the burn probability raster...")
terra::writeRaster(burn_prob, filename=paste0("./recalc_bp_", opt$foa_run, "_", opt$scenario, "_", opt$run_timepoint, ".tif"), overwrite=TRUE)

# Conditional FLPs: flame length probability given that it burned
cflps <- map(accum_flps, ~ .x / accum_bp)

# Unconditional FLPs: probability of this flame length occurring overall
flps <- map(cflps, ~ .x / num_seasons)

#Write cFLP and FLP rasters
flp_names <- names(categories)
for (i in seq_along(flp_names)) {
  name <- flp_names[i]
  log_message(paste0("Writing conditional flame length probability raster for flame lengths ", name))
  terra::writeRaster(cflps[[i]], filename = paste0("./recalc_cflp_", name, "_", opt$foa_run, "_", opt$scenario, "_", opt$run_timepoint, ".tif"), overwrite = TRUE)
  log_message(paste0("Writing unconditional flame length probability raster for flame lengths ", name))
  terra::writeRaster(flps[[i]], filename = paste0("./recalc_flp_", name, "_", opt$foa_run, "_", opt$scenario, "_", opt$run_timepoint, ".tif"), overwrite = TRUE)
}

# Clean up parallel backend
plan(sequential)
