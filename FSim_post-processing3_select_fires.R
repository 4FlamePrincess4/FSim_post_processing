library(tidyverse)
library(sp)
library(sf)
library(raster)

#Set the working directory to the specific outputs folder for the run
wd <- setwd("D:/WFSETP/FSim_Outputs/2022_baseline_t0/okwen_foa2c_r10/")

#######################################################################################
# NOTE: To run this code, you need to make sure the following FSim outputs are in the #
#       working directory: the FireSizeList.csv files, the All.tif files, the         #
#       ArrivalDays tifs, the FlameLengths tifs, and the Perimeters.sqlite files.     #
#######################################################################################

#STEP 1: Record run information below 
###############################################
foa_run <- "FOA2c_r10"
#scenario <- "LF2020"
run_timepoint <- "baseline_time0"
foa_lcp <- raster("../../../Data/Landfire2022_LCGs/FOA2c_Landfire2022/FOA2c_LCG_LF2022_FBFM40_230_120m.tif")
okwen_perimeter <- st_read("../../../Data/OkWen_shapefiles/FOA_shapefiles/OkWen_AllFOAs_60km_buffer/OkWen_cFOAs_Albers_60km_Buffer.shp")
foa_extent <- raster::extent(foa_lcp)
okwen_extent <- raster::extent(okwen_perimeter)
number_of_seasons <- 20000
#Use the below if you have an equal number of seasons for each part
number_of_parts <- 8
seasons_per_part <- c(rep(2500, number_of_parts))
#Use the below alternative if you have different numbers of seasons for each part
seasons_per_part <- c(5000, 7000, 2000, 1000, 5000)

#STEP 2: Combine fire lists from all four run parts then randomly select seasons
###############################################
#Read in run fire lists
firelist_files <- list.files(path=wd,
                             recursive=F,
                             pattern=".+FireSizeList.csv$",
                             full.names=T)

firelist_tables <- lapply(firelist_files, read.csv, header=TRUE)
firelists <- do.call(rbind, firelist_tables)
firelists <- firelists %>%
  mutate(Season_FireID = paste0(Season,"_",FireID))

#Solution to create part labels for parts of varying numbers of seasons
#First, make x number of vectors, where x is equal to the number of parts, 
# with a unique sequence of numbers corresponding to the seasons for each part
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
#Second, compare each fire's season to a list of seasons for each part and assign the 
# correct part label
#Initialize an empty vector and add it as a column to the firelists dataframe
Part <- vector("character",nrow(firelists))
firelists <- cbind(firelists,Part)
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

#Randomly identify seasons to be used for burn simulation in FFE-FVS
set.seed(42)
total_seasons <- 1:number_of_seasons
random_seasons <- sample(total_seasons, 10)
random_fires <- firelists %>% 
  dplyr::filter(Season %in% random_seasons)
random_fireIDs <- random_fires$FireID

#Write the random fires record to preserve all information associated with the selected seasons and fires
write_csv(random_fires, paste0("./", foa_run, "_selected_fires.csv"))

#STEP 2: Calculate summary statistics for each season and use these stats
# to select x number of seasons across the distribution of seasonal burned area
###############################################
#Specify the number of seasons to use for the daily emissions geodatabase
x <- 10
#Initialize a dataframe to hold the stats
season_stats <- data.frame(matrix(ncol=4,nrow=0,
                                  dimnames=list(NULL,
                                                c("Season",
                                                  "Acres_burned",
                                                  "Num_fires",
                                                  "Mu_fire_size"))))
#Calculate total area burned, total fires, and mean fire size for each season
for(each_season in unique(firelists$Season)){
  this_season <- dplyr::filter(firelists, Season == each_season)
  burned_area <- sum(this_season$Acres)
  tot_fires <- nrow(this_season)
  mu_firesize <- mean(this_season$Acres)
  season_stats[nrow(season_stats)+ 1,] = c(each_season,
                                           burned_area,
                                           tot_fires,
                                           mu_firesize)
}
write_csv(season_stats, paste0("./", foa_run, "_season_stats.csv"))

#Create a probability density function to grab percentile years
foa_areaburned_pdf <- density(season_stats$Acres_burned)
plot(foa_areaburned_pdf)
hist(season_stats$Acres_burned, breaks=40)
abline(v=areaburned_quantiles, col='red',lwd=1,lty='dashed')
#vector of quantiles
areaburned_quantiles <- stats::quantile(season_stats$Acres_burned, probs = seq(0, 1, (1/x)),names=FALSE,type=3)
print(areaburned_quantiles)
#Fetch the rows of season stats corresponding to the quantile values
quantile_season_stats <- season_stats[season_stats$Acres_burned %in% areaburned_quantiles,] 
print(quantile_season_stats)