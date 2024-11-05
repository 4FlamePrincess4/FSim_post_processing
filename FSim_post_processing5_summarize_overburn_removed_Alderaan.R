library(tidyverse)
library(raster)
library(terra)
library(tidyterra)
library(RSQLite)
library(furrr)
library(optparse)

#Set up input arguments with optparse
option_list = list(
  make_option(c("-w", "--working_directory"), type="character", default=NULL,
              help="working directory (mandatory)", metavar="character"),
  make_option(c("-f", "--season_fires_directory"), type="character", default=NULL,
              help="season fires directory (mandatory)", metavar="character"),
  make_option(c("-r", "--foa_run"), type="character", default=NULL,
              help="foa run label (mandatory)", metavar="character"),
  make_option(c("-s", "--scenario"), type="character", default=NULL,
              help="project scenario (mandatory)", metavar="character"),
  make_option(c("-t", "--run_timepoint"), type="character", default=NULL,
              help="timepoint for the scenario (mandatory)", metavar="character"),
  make_option(c("-n", "--number_of_seasons"), type="integer", default=NULL,
              help="total number of seasons (mandatory)", metavar="integer"),
  make_option(c("--seasons_in_part"), type="integer", default=NULL,
              help="number of seasons in a part", metavar="integer"),
  make_option(c("--number_of_parts"), type="integer", default=NULL,
              help="number of run parts", metavar="integer"),
  make_option(c("-j", "--seasons_per_part"), type="character", default=NULL,
              help="vector of number of seasons in a part", metavar="character")
)
# parse the command-line options
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#Set the working directory to the specific outputs folder for the run
setwd(opt$working_directory)
wd <- getwd()

#######################################################################################
# NOTE: To run this code, you need to make sure the following FSim outputs are in the #
#       working directory: the FireSizeList.csv files and the SeasonFires_merged_tifs #
#       directory.                                                                    #
#######################################################################################

#STEP 1: Record run information below 
###############################################
# calculate or parse seasons_per_part
if (is.null(opt$seasons_per_part)) {
  seasons_per_part <- rep(opt$seasons_in_part, opt$number_of_parts)
} else {
  seasons_per_part <- as.integer(unlist(strsplit(opt$seasons_per_part, ",")))
}
# Set the optparse variables as local variables to then pass to furr_options() for parallelization
foa_run <- opt$foa_run
scenario <- opt$scenario

#STEP 2: Combine fire lists from all four run parts
###############################################
#Read in run fire lists
firelist_files <- list.files(path=wd,
                             recursive=F,
                             pattern=".+FireSizeList.csv$",
                             full.names=T)
firelist_tables <- lapply(firelist_files, read.csv, header=TRUE)
#Combine the fire lists from all of the parts
firelists <- do.call(rbind, firelist_tables)
#Add a Season_FireID variable to distinguish fires with the same numbers but from different parts
firelists <- firelists %>%
  mutate(Season_FireID = paste0(Season,"_",FireID))

### Solution to create part labels for parts of varying numbers of seasons: ###

## First, make x number of vectors, where x is equal to the number of parts, 
## with a unique sequence of numbers corresponding to the seasons for each part.

#Grab the cumulative number of seasons per part
cumsum_seasons <- cumsum(seasons_per_part)
#remove the object "part_seasons" so that the ifelse statement in lines 56 through 61 works
remove(part_seasons)
part_seasons_list <- list()
for(j in 1:length(seasons_per_part)){
  if(exists("part_seasons")){
    part_seasons <- c((1+max(part_seasons)):(cumsum_seasons[j]))
    print(part_seasons)
  }else{
    part_seasons <- c(1:seasons_per_part[j])
  }
  part_seasons_li <- list(assign(paste0("pt",j,"_seasons_",min(part_seasons),"_",max(part_seasons)), part_seasons))
  part_seasons_list <- append(part_seasons_list,part_seasons_li)
}

## Second, compare each fire's season to a list of seasons for each part and assign the 
## correct part label.

#Initialize an empty vector and add it as a column to the firelists dataframe
Part <- vector("character",nrow(firelists))
Scenario <- rep(scenario, nrow(firelists))
firelists <- cbind(firelists,Part,Scenario)
#Sort the firelists dataframe by Season number
firelists <- firelists[order(firelists$Season),]
#For each fire in the firelists dataframe
for(fire_record in 1:nrow(firelists)){
  #And for each part in the list of part vectors (with Season number sequences)
  for(part_seasons_li in seq_along(part_seasons_list)){
    #Make the list item a vector again
    this_vec <- unlist(part_seasons_list[part_seasons_li])
    #Compare each fire Season number to the vector of seasons corresponding to the part
    # and, if the fire Season number is in the vector, add the part label to the Part column
    if(firelists$Season[fire_record] %in% this_vec){
      firelists$Part[fire_record] <- paste0(sprintf("pt%s",part_seasons_li))
    }
  }
}

write_csv(firelists, paste0("./", foa_run,"_",scenario, "_merged_firelists.csv"))

#STEP 3: Estimate area removed during the overburn process
#############################################################
tif_files <- list.files(opt$season_fires_directory, pattern = "\\.tif$", full.names=TRUE)

# Set up parallel backend for a Linux cluster
# Use `multicore` for single-node, or `cluster` with specified workers for multi-node
# Replace `<number_of_cores>` with the number of cores available
plan(multicore) # Or plan(cluster, workers = <number_of_cores>)

# Specify log file path
log_file <- paste0("./overburn_est_log_", foa_run, "_", scenario, ".log")

# Helper function to log messages to the log file
log_message <- function(message) {
  cat(message, "\n", file = log_file, append = TRUE)
}

# Define the processing function for a single tif
process_tif <- function(tif) {
  # Extract the season number from the filename using a regular expression
  season_number <- as.numeric(sub(".*Season([0-9]+)_.*", "\\1", basename(tif)))
  log_message(paste0("Processing season ", season_number, "..."))
  
  # Load the raster
  raster <- rast(tif, lyrs = 1)
  
  # Calculate number of fires after overburn
  num.fires.after.overburn <- length(na.omit(unique(values(raster))))
  
  # Filter the firelists data frame for the current season
  season_fires <- firelists %>%
    filter(Season == season_number)
  
  # Calculate IDs present in firelist but not in raster and vice versa
  firelist_IDs <- season_fires$FireID
  raster_IDs <- unique(na.omit(values(raster)))
  
  IDs_in_firelist_not_raster <- setdiff(firelist_IDs, raster_IDs)
  IDs_in_raster_not_firelist <- setdiff(raster_IDs, firelist_IDs)
  
  # Log information on IDs
  log_message(paste0("These fire IDs are in the firelist for season ", season_number, 
                     " but not in the merged season raster: ", IDs_in_firelist_not_raster))
  log_message(paste0("These fire IDs are in the merged season raster for season ", season_number, 
                     " but not in the season firelist: ", IDs_in_raster_not_firelist))
  
  # Check number of fires before and after overburn
  num.fires.before.overburn <- nrow(season_fires)
  
  if (num.fires.after.overburn == num.fires.before.overburn) {
    log_message(paste0("The number of fires after deleting overburn is the same (", num.fires.after.overburn, ")."))
  } else {
    log_message(paste0("The number of fires after deleting overburn is ", num.fires.after.overburn, 
                       " and the number of fires before deleting overburn is ", num.fires.before.overburn, "."))
  }
  
  # Count the number of non-NA pixels
  non_na_count <- sum(!is.na(values(raster)))
  log_message(paste0("A total of ", non_na_count, " pixels burned in season ", season_number, "."))
  
  # Return a data frame with the season and non-NA count
  data.frame(Season = season_number, area_burned_pixel_count = non_na_count)
}

# Apply the processing function in parallel to each tif file
pixels_burned <- future_map_dfr(tif_files, process_tif)

# Summarize area burned by season 
firelist_summary <- firelists %>%
  group_by(Season) %>%
  summarize(acres_burned = sum(Acres)) %>%
  left_join(pixels_burned, by = "Season") %>%
  mutate(acres_no_overburn = area_burned_pixel_count * 14400 / 4047,
         overburn_acres = pmax(acres_burned - acres_no_overburn, 0))

# Write the results to CSV
output_path <- paste0("./overburn_correction_by_season_", foa_run, "_", scenario, "_", run_timepoint, ".csv")
write_csv(firelist_summary, output_path)

# Clean up parallel backend
plan(sequential)
