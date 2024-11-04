library(tidyverse)
library(raster)
library(terra)
library(tidyterra)
library(RSQLite)

#Set the working directory to the specific outputs folder for the run
setwd("D:/WFSETP/FSim_Outputs/TTN_LF2016/okwen_foa1c_r16_LF2016/")
wd <- getwd()

#######################################################################################
# NOTE: To run this code, you need to make sure the following FSim outputs are in the #
#       working directory: the FireSizeList.csv files, the All.tif files, the         #
#       ArrivalDays tifs, the FlameLengths tifs, and the Perimeters.sqlite files.     #
#######################################################################################

#STEP 1: Record run information below 
###############################################
foa_run <- "FOA1c_r16"
scenario <- "LF2016"
run_timepoint <- "baseline_time0"
number_of_seasons <- 20000
#Use the below if you have an equal number of seasons for each part
number_of_parts <- 4
seasons_per_part <- c(rep(5000, number_of_parts))
#Use the below alternative if you have different numbers of seasons for each part
#seasons_per_part <- c(5000, 7000, 2000, 1000, 5000)

#STEP 2: Read in combined firelists
###############################################
firelists <- read_csv(paste0("./", foa_run,"_",scenario, "_merged_firelists.csv"))

#STEP 3: Estimate area removed during the overburn process
#############################################################
season_fires_folder <- "./SeasonFires_merged_tifs_LF2016/"

tif_files <- list.files(season_fires_folder, pattern = "\\.tif$", full.names=TRUE)

# Create an empty data frame to store the season number and non-NA pixel counts
pixels_burned <- data.frame(Season = numeric(), area_burned_pixel_count = numeric())

for(tif in tif_files){
  # Extract the season number from the filename using a regular expression
  # Assumes filenames follow the pattern "Season#_merged_IDs_ADs_FLs.tif"
  season_number <- as.numeric(sub(".*Season([0-9]+)_.*", "\\1", basename(tif)))
  print(paste0("Processing season ", season_number, "..."))
  # Load the raster
  raster <- rast(tif, lyrs=1)
  plot(raster)
  
  #print(na.omit(unique(values(raster))))
  num.fires.after.overburn <- length(na.omit(unique(values(raster))))
  
  season_fires <- firelists %>%
    filter(Season == season_number)
  
  firelist_IDs <- season_fires$FireID
  raster_IDs <- unique(na.omit(values(raster)))
  #print(paste0("Fire IDs in season ", season_number, " from the FireSizeList file:", ))
  #print(paste0("Fire IDs in season ", season_number, " from the ID raster layer:", ))
  
  IDs_in_firelist_not_raster <- setdiff(firelist_IDs,raster_IDs)
  print(paste0("These fire IDs are in the firelist for season ", season_number, 
               " but not in the merged season raster: ", IDs_in_firelist_not_raster))
  IDs_in_raster_not_firelist <- setdiff(raster_IDs, firelist_IDs)
  print(paste0("These fire IDs are in the merged season raster for season ", season_number, 
               " but not in the season firelist: ", IDs_in_raster_not_firelist))
  
  num.fires.before.overburn <- nrow(season_fires)
  
  if (num.fires.after.overburn == num.fires.before.overburn) {
    print(paste0("The number of fires after deleting overburn is the same (", num.fires.after.overburn, ")."))
  } else {
    print(paste0("The number of fires after deleting overburn is ", num.fires.after.overburn, 
                 " and the number of fires before deleting overburn is ", num.fires.before.overburn, "."))
  }
  
  #Count the number of non-NA pixels
  non_na_count <- sum(!is.na(values(raster)))
  print(paste0("A total of ", non_na_count, " pixels burned in season ", season_number, "."))
  
  # Append the season number and non-NA count to the data frame
  pixels_burned <- rbind(pixels_burned, data.frame(Season = season_number, area_burned_pixel_count = non_na_count))
}
view(pixels_burned)

#Summarize area burned by season 
firelist_summary <- data.frame(Season = numeric(), acres_burned = numeric())

for(season in firelists$Season){
  season_fires <- firelists %>%
    dplyr::filter(Season == season)
  sum_acres <- sum(season_fires$Acres)
  firelist_summary <- rbind(firelist_summary, data.frame(Season=season, acres_burned = sum_acres))
}
firelist_summary <- unique(firelist_summary)
view(firelist_summary)

#Append the information to the original firelist dataframe
firelist_summary <- dplyr::left_join(firelist_summary, pixels_burned)

#Subtract the overburn acres from OG acres as a new column
firelist_summary <- firelist_summary %>%
  mutate(acres_no_overburn = area_burned_pixel_count*14400/4047)%>%
  mutate(overburn_acres = acres_burned - acres_no_overburn)
firelist_summary$overburn_acres[firelist_summary$overburn_acres < 0] <- 0
view(firelist_summary)
write_csv(firelist_summary, paste0("./overburn_correction_by_season_", foa_run, "_", scenario, "_", run_timepoint, ".csv"))

#Summarize the overburn acres
summary(firelist_summary$overburn_acres)
hist(firelist_summary$overburn_acres, main="Acres of overburn removed each season", breaks = 50)
ggplot(firelist_summary, aes(x = "", y = overburn_acres)) +
  geom_violin(trim = FALSE, fill = "lightblue") +
  labs(title = "Acres of overburn removed each season", x = "Overburn", y = "Acres") +
  theme_minimal()+
  scale_y_continuous(limits=c(0,5))
num_greater_than_10acres <- sum(firelist_summary$overburn_acres > 10)
num_greater_than_50acres <- sum(firelist_summary$overburn_acres > 50)
num_greater_than_100acres <- sum(firelist_summary$overburn_acres > 100)
num_greater_than_500acres <- sum(firelist_summary$overburn_acres > 500)
num_greater_than_1000acres <- sum(firelist_summary$overburn_acres > 1000)
num_greater_than_10000acres <- sum(firelist_summary$overburn_acres > 10000)