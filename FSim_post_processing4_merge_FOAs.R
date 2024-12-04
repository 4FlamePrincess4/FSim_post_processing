library(tidyverse)
library(raster)
library(optparse)

# Set up input arguments with optparse
option_list = list(
  make_option(c("-w", "--working_directory"), type="character", default=NULL,
              help="working directory (mandatory)", metavar="character"),
  make_option(c("-s", "--scenario"), type="character", default=NULL,
              help="project scenario (mandatory)", metavar="character"),
  make_option(c("-t", "--run_timepoint"), type="character", default=NULL,
              help="timepoint for the scenario (mandatory)", metavar="character"),
  make_option(c("-e", "--study_area_lcp"), type="character", default=NULL,
              help="path for study area tif file to which all the FOA rasters should snap", metavar="character")
)

# Parse the command-line options
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Set the working directory
setwd(opt$working_directory)
wd <- getwd()

# Create a directory to hold the results
merged_dir <- paste0("./okwen_", opt$run_timepoint, "_", opt$scenario, "_merged_outputs")
dir.create(merged_dir, showWarnings = FALSE)

# Read the study area raster and extract CRS, origin, and extent
study_area <- raster(opt$study_area_lcp)
study_area_crs <- crs(study_area)
study_area_origin <- origin(study_area)
study_area_extent <- extent(study_area)

# Locate subdirectories containing the scenario in their names
target_dirs <- list.dirs(path = wd, recursive = FALSE, full.names = TRUE) %>%
  .[grepl(opt$scenario, .)]

if (length(target_dirs) == 0) {
  stop("No matching directories found for the specified scenario.")
}

# Patterns to match each burn probability raster
patterns <- c(
  ".*FullRun_bp\\.tif$",
  ".*FullRun_flp_0to2\\.tif$",
  ".*FullRun_flp_2to4\\.tif$",
  ".*FullRun_flp_4to6\\.tif$",
  ".*FullRun_flp_6to8\\.tif$",
  ".*FullRun_flp_8to12\\.tif$",
  ".*FullRun_flp_12plus\\.tif$"
)

# Process each pattern
for (pattern in patterns) {
  message(paste("Processing pattern:", pattern))
  
  # Collect all rasters matching the pattern from the target directories
  matched_files <- unlist(lapply(target_dirs, function(dir) {
    list.files(path = dir, pattern = pattern, full.names = TRUE)
  }))
  
  if (length(matched_files) == 0) {
    message(paste("No rasters found for pattern:", pattern))
    next
  }
  
  # Load and process the rasters
  processed_rasters <- lapply(matched_files, function(file) {
    r <- raster(file)
    # Match CRS, origin, and extent to the study area
    crs(r) <- study_area_crs
    origin(r) <- study_area_origin
    r <- raster::extend(r, study_area_extent, value = NA)
    return(r)
  })
  
  # Merge the processed rasters
  if (length(processed_rasters) > 1) {
    merged_raster <- do.call(raster::mosaic, c(processed_rasters, fun = sum, na.rm = TRUE, tolerance = 4))
  } else {
    merged_raster <- processed_rasters[[1]]
  }
  
  # Define the output filename
  output_filename <- paste0(
  merged_dir, "/okwen_", opt$run_timepoint, "_", opt$scenario, "_", 
  gsub("\\.\\*\\$", "", pattern), "_merged.tif"
)
  
  # Save the merged raster
  writeRaster(merged_raster, output_filename, format = "GTiff", overwrite = TRUE)
  message(paste("Merged raster saved to:", output_filename))
}
