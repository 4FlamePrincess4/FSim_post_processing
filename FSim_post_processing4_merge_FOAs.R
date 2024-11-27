library(tidyverse)
library(sf)
library(sp)
library(raster)
library(terra)
library(tidyterra)
library(optparse)

#Set up input arguments with optparse
option_list = list(
  make_option(c("-w", "--working_directory"), type="character", default=NULL,
              help="working directory (mandatory)", metavar="character"),
  make_option(c("-s", "--scenario"), type="character", default=NULL,
              help="project scenario (mandatory)", metavar="character"),
  make_option(c("-t", "--run_timepoint"), type="character", default=NULL,
              help="timepoint for the scenario (mandatory)", metavar="character"),
  make_option(c("-e","--study_area_lcp"), type="character", default=NULL,
              help="path for study area tif file to which all the FOA rasters should snap", metavar="character")
)
# parse the command-line options
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#Set the working directory to the specific outputs folder for the run
setwd(opt$working_directory)
wd <- getwd()

#######################################################################################
# NOTE: To run this code, you need to make sure the following FSim outputs are in the #
#       working directory: the FireSizeList.csv files, the All.tif files, the         #
#       ArrivalDays tifs, the FlameLengths tifs, and the Perimeters.sqlite files.     #
#######################################################################################
#Create a directory to hold the results 
merged_dir <- paste0("./okwen_", opt$run_timepoint, "_", opt$scenario, "_merged_outputs")
dir.create(merged_dir)

#Create a function that matches the crs, origin, and extent of the FOA rasters to that of the study area raster.
match_foas_to_study_area <- function(foa_raster_list, reference_crs, reference_origin, reference_extent){
  #Apply the operations to each raster in the list
  processed_rasters <- lapply(foa_raster_list, function(r){
    crs(r) <- reference_crs
    origin(r) <- reference_origin
    extended_r <- raster::extend(r, reference_extent, value = NA)
    return(extended_r)
    })
  return(processed_rasters)
  }

#Read in the study area template raster and save its crs and origin.
okwen_extent <- raster(opt$study_area_lcp)
okwen_crs <- crs(okwen_extent)
okwen_origin <- origin(okwen_extent)

#Locate subdirectories containing the scenario and timepoint in their names
target_dirs <- list.dirs(path= wd, recursive = FALSE, full.names = TRUE)

#Filter the subdirectories to only those that match the scenario and timepoint
target_dirs <- target_dirs[grep1(opt$scenario, target_dirs)]
#Use the below if you're working with a timepoint run for the ongoing project
#target_dirs <- target_dirs[grep1(opt$scenario, target_dirs) & grep1(opt$run_timepoint, target_dirs)]

#Check if there are matching directories
if(length(target_dirs) == 0){
  stop("No matching directories found for the specified scenario and timepoint.")
  }

#Patterns to match each burn probability raster
patterns <- c(
  ".*FullRun_bp\\.tif$",
  ".*FullRun_flp_0to2\\.tif$",
  ".*FullRun_flp_2to4\\.tif$",
  ".*FullRun_flp_4to6\\.tif$",
  ".*FullRun_flp_6to8\\.tif$",
  ".*FullRun_flp_8to12\\.tif$",
  ".*FullRun_flp_12plus\\.tif$"
)

#Define a function to search within restricted subdirectories
function(pattern) {
  list.files(
    path = target_dir, # Search within the restricted subdirectory
    pattern = pattern,
    full.names = TRUE
  )
}

# Search for files matching the patterns in the target directories
matched_files <- lapply(target_dirs, function(dir) {
  lapply(patterns, function(pattern) {
    list.files(path = dir, pattern = pattern, full.names = TRUE)
  })
})
# Flatten the list of matched files into a single vector
matched_files <- unlist(matched_files)

# Print matched files for verification
print(matched_files)
  
# Convert matched_files (file paths) into a list of raster objects
foa_raster_list <- lapply(matched_files, raster)
  
# Match CRS, origin, and extent to the study area raster
processed_rasters <- match_foas_to_study_area(foa_raster_list, okwen_crs, okwen_origin, okwen_extent)
  
# Prepare arguments for raster::mosaic
names(processed_rasters)[1:2] <- c('x', 'y') # Rename for mosaic
processed_rasters$fun <- sum
processed_rasters$na.rm <- TRUE
processed_rasters$tolerance <- 4
  
# Use raster::mosaic to merge the rasters
merged_raster <- do.call(raster::mosaic, processed_rasters)
  
# Define the output filename based on the pattern
output_filename <- paste0(
  merged_dir, "/okwen_", opt$run_timepoint, "_", opt$scenario, "_", 
  gsub(".+FullRun_", "", gsub("\\.tif$", "", pattern)), "_merged.tif"
)
  
# Save the merged raster
writeRaster(merged_raster, output_filename, format = "GTiff", overwrite = TRUE)
message(paste("Merged raster saved to:", output_filename))
