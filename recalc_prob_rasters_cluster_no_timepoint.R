library(tidyverse)
library(terra)
library(tidyterra)
library(optparse)

#Set up input arguments with optparse
option_list = list(
  make_option(c("-w", "--working_directory"), type="character", default=NULL,
              help="working directory (mandatory)", metavar="character"),
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

#List SeasonFire raster paths
season_fire_files <- list.files(
  path = paste0("./SeasonFires_merged_tifs_", opt$scenario),
  pattern = ".tif$", full.names=TRUE
)

# Use a virtual raster to avoid loading all into memory
season_fire_rasters <- vrt(season_fire_files)

#Compute burn probability
prop_non_na <- app(season_fire_rasters, function(x) mean(!is.na(x)), filename=paste0("./recalc_", opt$foa_run, "_", opt$scenario, "_", opt$run_timepoint, "_burn_probability.tif"), overwrite=TRUE)
names(prop_non_na) <- "burn_probability"

# Create a temporary directory for storing intermediate tifs
temp_dir <- file.path(wd, "temp_tifs")
dir.create(temp_dir)

# Process flame length probabilities efficiently
categories <- list(
  "0_2"   = c(0, 2),
  "2_4"   = c(2, 4),
  "4_6"   = c(4, 6),
  "6_8"   = c(6, 8),
  "8_12"  = c(8, 12),
  "12plus" = c(12, Inf)
)

# Process flame length band (3rd band)
for (cat in names(categories)) {
  bounds <- categories[[cat]]
  
  # Process each raster file one at a time
  category_count_list <- list()
  
  for (i in seq_along(season_fire_files)) {
    r <- rast(season_fire_files[i], lyr=3) # Load only 3rd band
    
    # Fix: Use a function that returns a single value per pixel
    category_count_raster <- app(r, function(x) { 
      ifelse(is.na(x), NA, sum(x >= bounds[1] & x < bounds[2])) 
    })
    
    # Save the result and add to the list
    filename_temp <- file.path(temp_dir, paste0("temp_", i, ".tif"))
    writeRaster(category_count_raster, filename_temp, overwrite=TRUE)
    category_count_list[[i]] <- rast(filename_temp)
  }
  
  # Merge counts across all files using tapp()
  category_count_raster <- tapp(rast(category_count_list), fun=sum, 
                                filename=paste0("recalc_", opt$foa_run, "_", opt$scenario, "_", opt$run_timepoint, "_cflp_flame_length_", cat, ".tif"),
                                overwrite=TRUE)
  
  # Compute probability
  probability_raster <- category_count_raster / prop_non_na
  writeRaster(probability_raster, paste0("recalc_", opt$foa_run, "_", opt$scenario, "_", opt$run_timepoint, "_flp_flame_length_", cat, ".tif"), overwrite=TRUE)
  
  # Cleanup: Delete all temporary files and remove the directory
  unlink(temp_dir, recursive=TRUE)
}
