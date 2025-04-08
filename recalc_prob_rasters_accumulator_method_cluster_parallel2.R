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

calc_prob_w_accumulator_safe <- function(season_fire_path, category_values, foa_lcp_path, bp_template_path, fl_template_paths, output_dir) {
  tryCatch({
    #Load rasters inside of worker
  foa_lcp <- terra::rast(opt$foa_lcp_path, lyrs=1)
  foa_lcp <- terra::unwrap(foa_lcp)

  #Create an empty bp accumulator raster
  accum_bp <- terra::rast(foa_lcp)
  terra::values(accum_bp) <- NA

  #Create empty accumulator rasters for the flame length probability rasters
  fl_accumulators <- lapply(names(categories), function(cat) {
    r <- setValues(foa_lcp, NA)
    names(r) <- paste0("acc_", cat)
    return(r)
  })
  names(fl_accumulators) <- names(categories)

    # Run original function (assuming it returns SpatRasters)
    results <- calc_prob_w_accumulator(
      season_fire_path = season_fire_path,
      categories = category_values,
      foa_lcp = foa_lcp,
      accum_bp_template = accum_bp_template,
      fl_accumulators_template = fl_accumulators_template
    )

    # Write results to disk (return paths instead of SpatRaster)
    season_id <- tools::file_path_sans_ext(basename(season_fire_path))
    bp_path <- file.path(output_dir, paste0("bp_accum_", season_id, ".tif"))
    fl_paths <- purrr::imap(results$fl_accums, ~ {
      file.path(output_dir, paste0("fl_", .y, "_accum_", season_id, ".tif"))
    })

    terra::writeRaster(results$accum_bp, bp_path, overwrite = TRUE)
    purrr::walk2(results$fl_accums, fl_paths, ~ terra::writeRaster(.x, .y, overwrite = TRUE))

    return(list(bp_path = bp_path, fl_paths = fl_paths))
  }, error = function(e) {
    msg <- paste(Sys.time(), "- Error in", season_fire_path, ":", e$message)
    cat(msg, "\n", file = "error_log.txt", append = TRUE)
    return(NULL)
  })
}


calc_prob_w_accumulator <- function(season_fire_path, categories, fl_accumulators, output_dir){
   
  #Read in the SeasonFire raster (band 3 - flame lengths)
  seasonfire_FLs <- terra::rast(season_fire_path, lyr=3)
  #Set CRS to match foa_lcp
  terra::crs(seasonfire_FLs) <- terra::crs(foa_lcp)
  #Extend raster to match extent of foa_lcp
  seasonfire_FLs <- terra::extend(seasonfire_FLs, terra::ext(foa_lcp), snap="near")
  
  #Create a mask of all pixels that burned at any flame length
  burned_mask <- !is.na(seasonfire_FLs)
  #Add the mask to the bp accumulator raster
  accum_bp <- accum_bp + burned_mask
  
  #Convert the flame length raster to integers for binning into flame length cats
  seasonfire_FLs_int <- round(seasonfire_FLs)
  
  #Generate one mask for each category
  for (cat in names(categories)) {
    bounds <- categories[[cat]]
    mask <- seasonfire_FLs_int >= bounds[1] & seasonfire_FLs_int < bounds[2]
    
    # Update accumulator: where mask is TRUE, add 1; where NA, leave as is
    fl_accumulators[[cat]] <- ifel(mask, 
                                   ifel(!is.na(fl_accumulators[[cat]]), fl_accumulators[[cat]] + 1, 1), 
                                   fl_accumulators[[cat]])
  }
  #Add the burn probability raster to the front of the list
  prob_rast_list <- c(list(accum_bp = accum_bp), fl_accumulators)
  #Return the full list
  return(prob_rast_list)
}

#List the SeasonFire raster paths
season_fire_files <- list.files(path = paste0(wd,"/SeasonFires_merged_tifs_", opt$scenario),
                                pattern = ".tif$", full.names=TRUE)
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

#Read in each seasonfire raster and update the accumulators
# Remember arguments must be in this order: season_fire_path, categories, fl_accumulators
log_message("Calculating accumulator bp and flp rasters...")

# Output directory (safe to share across workers if filenames are unique)
output_dir <- file.path(tempdir(), "accumulator_outputs")
dir.create(output_dir, showWarnings = FALSE)

for(season in 1:length(season_fire_files)){
  season_fire_path <- season_fire_files[season]
  result <- calc_prob_w_accumulator(season_fire_path, categories, fl_accumulators, output_dir)
}

#Separate the list of results into discrete rasters
accum_bp <- result$accum_bp
accum_flp1 <- result$acc_0_2
accum_flp2 <- result$acc_2_4
accum_flp3 <- result$acc_4_6
accum_flp4 <- result$acc_6_8
accum_flp5 <- result$acc_8_12
accum_flp6 <- result$acc_12plus
rm(result)

#Rename each of the rasters
names(accum_bp) <- "recalc_bp"
names(accum_flp1) <- "recalc_flp_0to2ft"
names(accum_flp2) <- "recalc_flp_2to4ft"
names(accum_flp3) <- "recalc_flp_4to6ft"
names(accum_flp4) <- "recalc_flp_6to8ft"
names(accum_flp5) <- "recalc_flp_8to12ft"
names(accum_flp6) <- "recalc_flp_12plusft"

#Calculate the burn probability raster
burn_prob <- accum_bp/num_seasons
#Write the burn probability raster
log_message("Writing the burn probability raster...")
writeRaster(burn_prob, filename=paste0("./recalc_bp_", foa_run, "_", scenario, ".tif"), overwrite=TRUE)
rm(burn_prob)

#Calculate the conditional flp rasters
cflp1 <- accum_flp1/accum_bp
cflp2 <- accum_flp2/accum_bp
cflp3 <- accum_flp3/accum_bp
cflp4 <- accum_flp4/accum_bp
cflp5 <- accum_flp5/accum_bp
cflp6 <- accum_flp6/accum_bp
rm(accum_bp, accum_flp1, accum_flp2, accum_flp3, accum_flp4, accum_flp5, accum_flp6)

#Write the conditional flp rasters
log_message(paste0("Writing conditional flame length probability rasters ..."))
writeRaster(cflp1, filename=paste0("./recalc_cflp_0to2ft_", foa_run, "_", scenario, ".tif"), overwrite = TRUE)
writeRaster(cflp2, filename=paste0("./recalc_cflp_2to4ft_", foa_run, "_", scenario, ".tif"), overwrite = TRUE)
writeRaster(cflp3, filename=paste0("./recalc_cflp_4to6ft_", foa_run, "_", scenario,".tif"), overwrite = TRUE)
writeRaster(cflp4, filename=paste0("./recalc_cflp_6to8ft_", foa_run, "_", scenario, ".tif"), overwrite = TRUE)
writeRaster(cflp5, filename=paste0("./recalc_cflp_8to12ft_", foa_run, "_", scenario,".tif"), overwrite = TRUE)
writeRaster(cflp6, filename=paste0("./recalc_cflp_12plusft_", foa_run, "_", scenario, ".tif"), overwrite = TRUE)

#Calculate the unconditional flp rasters
flp1 <- cflp1/num_seasons
flp2 <- cflp2/num_seasons
flp3 <- cflp3/num_seasons
flp4 <- cflp4/num_seasons
flp5 <- cflp5/num_seasons
flp6 <- cflp6/num_seasons

#Write the unconditional flp rasters
log_message(paste0("Writing unconditional flame length probability rasters ..."))
writeRaster(flp1, filename=paste0("./recalc_flp_0to2ft_", foa_run, "_", scenario, ".tif"), overwrite = TRUE)
writeRaster(flp2, filename=paste0("./recalc_flp_2to4ft_", foa_run, "_", scenario, ".tif"), overwrite = TRUE)
writeRaster(flp3, filename=paste0("./recalc_flp_4to6ft_", foa_run, "_", scenario, ".tif"), overwrite = TRUE)
writeRaster(flp4, filename=paste0("./recalc_flp_6to8ft_", foa_run, "_", scenario, ".tif"), overwrite = TRUE)
writeRaster(flp5, filename=paste0("./recalc_flp_8to12ft_", foa_run, "_", scenario, ".tif"), overwrite = TRUE)
writeRaster(flp6, filename=paste0("./recalc_flp_12plusft_", foa_run, "_", scenario, ".tif"), overwrite = TRUE)
rm(flp1, flp2, flp3, flp4, flp5, flp6, cflp1, cflp2, cflp3, cflp4, cflp5, cflp6)
gc()
