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
# parse the command-line options
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#Set the working directory to the specific outputs folder for the run
setwd(opt$working_directory)
wd <- getwd()

#List the SeasonFire raster paths
season_fire_files <- list.files(path = paste0("./SeasonFires_merged_tifs_", opt$scenario, "_", opt$run_timepoint),
                                pattern = ".tif$", full.names=TRUE)

#Create a SpatRaster dataset
season_fire_rasters <- rast(lapply(season_fire_files, function(f) rast(f, lyr=1)))

#Compute the propoprtion of non-NA values per cell to get at burn probability (# times burned/total)
prop_non_na <- app(season_fire_rasters, function(x) mean(!is.na(x)))
names(prop_non_na) <- "burn_probability"

#Write this burn probability raster
writeRaster(prop_non_na, paste0("./recalc_", opt$foa_run, "_", opt$scenario, "_", opt$run_timepoint, "_burn_probability.tif"))

#Let's do the same for the flame length probabilities
# Define flame length categories
categories <- list(
  "0_2"   = c(0, 2),
  "2_4"   = c(2, 4),
  "4_6"   = c(4, 6),
  "6_8"   = c(6, 8),
  "8_12"  = c(8, 12),
  "12plus" = c(12, Inf)
)

# Extract the third band (flame length data)
flame_length_rasters <- rast(lapply(season_fire_files, function(f) rast(f, lyr = 3)))

# Compute the total number of times a cell is not NA
non_na_count <- app(flame_length_rasters, function(x) sum(!is.na(x)))

# Compute probability rasters for each category
probability_rasters <- lapply(names(categories), function(cat) {
  bounds <- categories[[cat]]
  
  # Count how many times each cell falls into this category
  category_count_raster <- app(flame_length_rasters, function(x) {
    sum(x >= bounds[1] & x < bounds[2], na.rm = TRUE)
  })
  
  # Compute proportion raster
  probability_raster <- category_count_raster / non_na_count
  
  return(probability_raster)
})

# Calculate the unconditional flp rasters and save both conditional & unconditional versions
for(i in seq_along(probability_rasters)) {
  writeRaster(probability_rasters[[i]], paste0("recalc_", foa_run, "_", scenario, "_", opt$run_timepoint, "_cflp_flame_length_", names(categories)[i], ".tif"), overwrite = TRUE)
  unconditional_flp <- probability_rasters[[i]]*prop_non_na
  writeRaster(unconditional_flp, paste0("recalc_", foa_run, "_", scenario, "_", opt$run_timepoint, "_flp_flame_length_", names(categories)[i], ".tif"), overwrite = TRUE)
}
