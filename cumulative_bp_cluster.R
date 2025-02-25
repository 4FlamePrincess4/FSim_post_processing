library(tidyverse)
library(terra)
library(tidyterra)
library(optparse)

#Set up input arguments with optparse
option_list = list(
  make_option(c("-i", "--foa_lcp_path"), type="character", default=NULL,
              help="foa lcp path (mandatory)", metavar="character"),
  make_option(c("-w", "--working_directory"), type="character", default=NULL,
              help="working directory (mandatory)", metavar="character"),
  make_option(c("-r", "--foa_run"), type="character", default=NULL,
              help="foa run label (mandatory)", metavar="character"),
  make_option(c("-s", "--scenario"), type="character", default=NULL,
              help="project scenario (mandatory)", metavar="character"),
  make_option(c("-t", "--run_timepoint"), type="character", default=NULL,
              help="timepoint for the scenario (mandatory)", metavar="character"),
  make_option(c("-f", "--first_season"), type="integer", default=NULL,
              help="first season (mandatory)", metavar="integer"),
  make_option(c("-l", "--last_season"), type="integer", default=NULL,
              help="last season (mandatory)", metavar="integer"),
  make_option(c("-p", "--num_random_pixels"), type="integer", default=NULL,
              help="number of random pixels to sample (mandatory)", metavar="integer")
)
# parse the command-line options
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#Set the working directory to the specific outputs folder for the run
setwd(opt$working_directory)
wd <- getwd()

set.seed(42)

#Read in the first layer of the LCP
lcp <- rast(opt$foa_lcp_path, lyr=1)
#Replace LCP values with 0 
terra::values(lcp) <- 0

#Generate a sample of 100 random pixels from the total pixels in the lcp
pixel_count <- ncell(lcp)
random_pixels <- sample(1:pixel_count, opt$num_random_pixels)

#Initialize a dataframe to hold the results
cumulative_bp_df <- data.frame(matrix(numeric(0), nrow=0, ncol=opt$num_random_pixels))
colnames(cumulative_bp_df) <- paste0("pix", 1:opt$num_random_pixels)

#Create a for loop that iterates over the SeasonFire rasters and cumulatively calculates the bp
seasons <- opt$first_season:opt$last_season
for(season in seasons){
  #Read in SeasonFire raster
  this_season_rast <- rast(paste0(wd, "/SeasonFires_merged_tifs_", opt$scenario, "/Season", season, "_merged_IDs_ADs_FLs.tif"), lyr=1)
  # Create a binary mask where overlapping pixels are 1, else 0
  binary_mask <- ifel(!is.na(this_season_rast), 1, 0)
  # Add the binary mask to lcp, ensuring only one increment per overlapping pixel
  lcp <- lcp + binary_mask
  #Grab the pixel values corresponding to the random number vector
  random_raw_values <- terra::values(lcp)[random_pixels]
  #Divide by the current iteration number
  random_bp_values <- random_raw_values/season
  #Append the burn probabilities to the dataframe
  cumulative_bp_df[nrow(cumulative_bp_df) + 1,] <- as.numeric(random_bp_values)
}

#Add an iteration column to the dataframe
cumulative_bp_df$iteration <- 1:nrow(cumulative_bp_df)

#Reshape the dataframe to long format for ggplot
long_df <- cumulative_bp_df %>%
  pivot_longer(cols = starts_with("pix"),
               names_to = "pixel",
               values_to = "burn_prob")

#Export the dataframe
write_csv(long_df, paste0("./cumulative_bp_", opt$num_random_pixels, "_pixels_", opt$foa_run, "_", opt$scenario, ".csv"))

