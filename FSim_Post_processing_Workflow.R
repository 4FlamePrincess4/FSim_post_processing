library(tidyverse)
library(sf)
library(sp)
library(raster)
library(terra)
library(tidyterra)
library(RSQLite)

#Set the working directory to the specific outputs folder for the run
wd <- setwd("E:/Wildland_Fire_Smoke_Emissions_Tradeoff_Project/FSim_Outputs/baseline_t0/okwen_foa2c_r10/")

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
                             recursive=T,
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

#STEP 3: Merge (average) probability rasters from all parts of each run
###############################################
#Read in the "All" tifs, which are multiband rasters of outputs that include
## 1) burn probability, 2) conditional flame length probability (0-2 ft), 3) conditional flame
## length probability (2-4 ft), 4) CFL prob (4-6 ft), 5) CFL prob (6-8 ft), 6) CFL prob (8-12 ft),
## 7) CFL prob (12+ ft), 8) Conditional Flame Length (ft). 
all_tifs_files <- list.files(path=wd,
                             recursive=T,
                             pattern=".+All.tif$",
                             full.names=T)
alltifs_stacks <- lapply(all_tifs_files, stack)

#Create a data frame identifying each raster in the All stack and the corresponding band number.
layers <- data.frame(layer_name=c("bp","flp_0to2","flp_2to4",
                                  "flp_4to6","flp_6to8","flp_8to12",
                                  "flp_12plus"),layer_num=c(seq(from=1,to=7)))

#Initialize lists to hold tifs in the for loop
p_tifs <- list()
#For each of the rasters in a stack
for(this_layer in 1:nrow(layers)){
  #For each stack in the list of stacks
  for(this_stack in seq_along(alltifs_stacks)){
    #Isolate the raster layer of interest
    this_raster <- raster::subset(alltifs_stacks[[this_stack]], layers$layer_num[this_layer])
    #Also isolate the burn probability raster to make the conditional flp rasters unconditional
    bp_raster <- raster::subset(alltifs_stacks[[this_stack]], 1)
    #If the raster layer is one of the cflp layers, multiply by the bp raster
    uncond_p_raster <- if(between(layers$layer_num[this_layer],2,7)){
      this_raster*bp_raster
    }else{
      this_raster
    }
    #plot(uncond_p_raster)+title(paste0("unconditional ",layers$layer_name[this_layer],this_stack))
    #Print indicator
    print(paste0("Processing ",layers$layer_name[this_layer]," of part ",
                 this_stack,"..."))
    #Assign the raster to a variable with a name unique to the raster layer and stack,
    # and make it a list item.
    p_tif <- list(assign(paste0(layers$layer_name[this_layer],"_",this_stack), uncond_p_raster))
    #Append the list item to a list of rasters
    p_tifs <- append(p_tifs, p_tif)
  }
  #Take the average of all the rasters in the list
  ptif_stack <- do.call(stack,p_tifs)
  print("Processing the weighted mean of the stack of part rasters...")
  mean_p <- raster::weighted.mean(ptif_stack,seasons_per_part,na.rm=TRUE)
  #Plot the averaged raster
  plot(mean_p)+title(paste0(foa_run,"_",run_timepoint,"_FullRun_",layers$layer_name[this_layer]))
  #Write the averaged raster to the working directory
  writeRaster(mean_p,paste0(foa_run,"_",run_timepoint,"_FullRun_",layers$layer_name[this_layer]),
               format="GTiff", overwrite=TRUE)
  #Re-initialize the p_tifs list
  p_tifs <- list()
}

#STEP 4: Read all flame length tifs, filter by the randomly selected seasons & export
########################################################################################################################
# NOTE: Before you can run this code, you need to give the fire tif files unique names indicating their run and part.  #
#   To do this, run the PowerShell script "rename_tifs.ps1", which you can download from the FSim_scripts Google Drive #
#   folder. You need to place the script in the folder, then you can right click and select the option to run the      #
#   script in PowerShell. The script will be renamed as well, but this doesn't matter. I suggest copying the original  #
#   to each folder and running the script before moving on to pasting into the next folder. The script is a small file.#
########################################################################################################################
#Read in flame length tifs
# First, list the file names
FL_tif_list <- list.files(path=wd,
                          recursive=T,
                          pattern=".+FlameLengths.+.tif$",
                          full.names=T)
# Next, select the flame length tifs selected for burn simulation in FFE-FVS
#Initialize a list to hold labels
raster_labels <- list()
#Remove a flame length raster stack from a previous run, if applicable
#This is necessary because of the if-else statement in lines 145-151.
remove(selected_seasons_fl_raster_stack)
#For each season
for(each_season in unique(random_fires$Season)){
  #Subset the random_fires dataframe to the season in question
  this_season_fires <- random_fires %>%
    filter(Season == each_season)
  #Grab the fire IDs vector and make sure they are strings
  this_season_fireIDs <- this_season_fires$FireID
  this_season_fireIDs <- as.character(this_season_fireIDs)
  #Grab the part and make sure they are strings
  this_season_pt <- this_season_fires$Part
  this_season_pt <- as.character(this_season_pt)
  #Make a vector with the foa_run label of the same length as the fire IDs and part
  this_season_foa_run <- rep(foa_run, length(this_season_fireIDs))
  #Make a vector of file names with which to search the list of all Flame Length tifs
  this_season_fires_filenames <- paste0(this_season_foa_run, "_",
                                        this_season_pt, "_FlameLengths_FireID_",
                                        this_season_fireIDs, ".tif")
  #Initialize a new empty list to hold the tifs
  this_season_fires_fl_db <- list()
  #For each flame length tif file name
  for(each_fire in seq_along(this_season_fires_filenames)){
    #Search the list of all Flame Length tifs and grab the full path
    this_fire_filename <- grep(this_season_fires_filenames[each_fire],
                               FL_tif_list, value = TRUE, fixed=TRUE)
    #Print indicator for progress
    print(this_fire_filename)
    #Read in the flame length tif
    this_fire_tif <- raster::raster(this_fire_filename)
    #Plotting takes a long time but you can do this if you want another indicator
    #plot(this_fire_tif)
    #Add the fire tif to the list of tifs
    this_season_fires_fl_db <- append(this_season_fires_fl_db, this_fire_tif)

  }
  #Add the names and arguments needed for the raster::mosaic function.
  #We want to make one raster with all the fire flame length tifs, so we use the sum function.
  #Note: if a season happens to have overlapping fires, you'll run into problems here.
  #Very occasionally a season has only one fire, which throws an error with the mosaic function.
  # Therefore, I added an ifelse statement.
  if(length(this_season_fires_fl_db) > 1){
    names(this_season_fires_fl_db)[1:2] <- c('x','y')
    this_season_fires_fl_db$fun <- sum
    this_season_fires_fl_db$na.rm <- TRUE
    #Add the flame length rasters together to make one raster for each season
    this_season_fires_fl_raster <- do.call(raster::mosaic, this_season_fires_fl_db)
    #Plot the resulting raster
    plot(this_season_fires_fl_raster) + title(paste0(each_season))
  }else{
    #Plot the resulting raster
    this_season_fires_fl_raster <- this_season_fires_fl_db[[1]]
    plot(this_season_fires_fl_raster) + title(paste0(each_season))
  }
  #The rasters all have to be the same extent in order to stack properly
  this_season_fires_fl_raster <- raster::extend(this_season_fires_fl_raster, 
                                                okwen_extent, value=NA)
  this_season_fires_fl_raster <- raster::crop(this_season_fires_fl_raster, 
                                              foa_extent)
  #Confirm the result if you want
  #plot(this_season_fires_fl_raster)
  #If you don't already have a raster stack, make it. Then add your raster layers.
  if(exists("selected_seasons_fl_raster_stack")){
    selected_seasons_fl_raster_stack <- addLayer(selected_seasons_fl_raster_stack, 
                                                 this_season_fires_fl_raster)
  }else{
    selected_seasons_fl_raster_stack <- this_season_fires_fl_raster
  }
  #Label the rasters so you know which season is which once the rasters are stacked.
  raster_labels <- append(raster_labels, each_season)
}
#Now, write each raster in the stack and label each with the season selected.
for(season in seq_along(selected_seasons_fl_raster_stack)){
  this_season_raster <- subset(selected_seasons_fl_raster_stack, 
                               season, drop=TRUE)
  print(this_season_raster)
  writeRaster(this_season_raster,
              paste0("./Season", raster_labels[season],"_Fires_FlameLength.tif"), 
              format = "GTiff", overwrite=TRUE)
}

#STEP 5: Select and combine the ArrivalDay tifs for selected years and fires
############################################################################
# First, list the file names
AD_tif_list <- list.files(path=wd,
                          recursive=T,
                          pattern=".+ArrivalDays.+.tif$",
                          full.names=T,
                          all.files=T)
# The comments for ArrivalDays tif processing mirror those for the FlameLength tif code in Step 3
raster_labels <- list()
remove(selected_seasons_ad_raster_stack)
for(each_season in unique(random_fires$Season)){
  this_season_fires <- random_fires %>%
    filter(Season == each_season)
  this_season_fireIDs <- this_season_fires$FireID
  this_season_fireIDs <- as.character(this_season_fireIDs)
  this_season_pt <- this_season_fires$Part
  this_season_pt <- as.character(this_season_pt)
  this_season_foa_run <- rep(foa_run, length(this_season_fireIDs))
  this_season_fires_filenames <- paste0(this_season_foa_run, "_", this_season_pt, 
                                        "_ArrivalDays_FireID_",
                                        this_season_fireIDs, ".tif")
  this_season_fires_ad_db <- list()
  for(each_fire in seq_along(this_season_fires_filenames)){
    this_fire_filename <- grep(this_season_fires_filenames[each_fire],
                               AD_tif_list, value = TRUE, fixed=TRUE)
    this_fire_tif <- raster::raster(this_fire_filename)
    #plot(this_fire_tif)
    this_season_fires_ad_db <- append(this_season_fires_ad_db, this_fire_tif)
    
  }
  if(length(this_season_fires_ad_db) > 1){
    names(this_season_fires_ad_db)[1:2] <- c('x','y')
    this_season_fires_ad_db$fun <- sum
    this_season_fires_ad_db$na.rm <- TRUE
    this_season_fires_ad_raster <- do.call(raster::mosaic, this_season_fires_ad_db)
    plot(this_season_fires_ad_raster)
    title(paste0(each_season))
  }else{
    this_season_fires_ad_raster <- this_season_fires_ad_db[[1]]
    plot(this_season_fires_ad_raster)
    title(paste0(each_season))
  }
  this_season_fires_ad_raster <- raster::extend(this_season_fires_ad_raster,
                                                okwen_extent, value=NA)
  this_season_fires_ad_raster <- raster::crop(this_season_fires_ad_raster,
                                              foa_extent)
  plot(this_season_fires_ad_raster) + title(paste0(each_season))
  if(exists("selected_seasons_ad_raster_stack")){
    selected_seasons_ad_raster_stack <- addLayer(selected_seasons_ad_raster_stack,
                                                 this_season_fires_ad_raster)
  }else{
    selected_seasons_ad_raster_stack <- this_season_fires_ad_raster
  }
  raster_labels <- append(raster_labels, each_season)
}

for(season in seq_along(selected_seasons_ad_raster_stack)){
  this_season_raster <- subset(selected_seasons_ad_raster_stack, 
                               season, drop=TRUE, overwrite=TRUE)
  print(this_season_raster)
  writeRaster(this_season_raster,
              paste0("./Season", raster_labels[season],"_ArrivalDays.tif"), 
              format = "GTiff", overwrite=TRUE)
}

#STEP 6: Export perimeters for all fires
#############################################################
#Read in fire perimeters
perimeter_dbs <- list.files(path=wd,
                        recursive=T,
                        pattern="_Perimeters.sqlite",
                        full.names=T)
# #This is a general solution to get the FOA name for naming the file
# foa_run <- gsub(".okwen_foa3c_r3/*", "", gsub("*_Perimeters.sqlite", "", perimeter_dbs[1]))

i=1

#Loop through each perimeter, manipulate, and then merge all into one file

for (i in seq_along(perimeter_dbs)){
  
  #Work with the sqlite data ----
  
  #connect to the database
  con <- dbConnect(drv=RSQLite::SQLite(), dbname=perimeter_dbs[i])
  
  #list all tables, besides the unnecessary "sqlite_sequence"
  tables <- dbListTables(con)[ dbListTables(con) !="sqlite_sequence"]
  
  #Initiate the sql_dataframes
  sql_dataframes <- vector("list", length=length(tables))
  
  ## create a dataframe for each table
  for (j in seq_along(tables)) {
    assign(paste0("df",j), dbGetQuery(conn=con, 
                                      statement=paste("SELECT * FROM '", 
                                                      tables[[j]], "'", sep="")))
  }
  
  #create a shapefile from the dataframes ----

  
  df2$GEOMETRY <- sf::st_as_sfc(structure(as.list(df2$GEOMETRY),class="blob"),
                                        crs = df3$srtext)
  
  df2 <- st_as_sf(df2)

  #save the file
  # st_write(df2,
  #             paste0("./", "perimeters_", foa_run),
  #             driver= "ESRI Shapefile", append=TRUE)
  st_write(df2,
           paste0("./", "perimeters_", foa_run,".geojson"),
           append=TRUE)
  # st_write(df2,
  #          paste0("./", "perimeters_", run_timepoint),
  #          driver= "ESRI Shapefile", append=TRUE)
}
dbDisconnect()

#STEP 7: Export ignitions for all fires
#############################################################
#Read in fire ignitions
point_dbs <- list.files(path=wd,
                            recursive=T,
                            pattern="_Ignitions.sqlite",
                            full.names=T)
# #This is a general solution to get the FOA name for naming the file
# foa_run <- gsub(".okwen_foa3c_r3/*", "", gsub("*_Perimeters.sqlite", "", perimeter_dbs[1]))

i=1

#Loop through each perimeter, manipulate, and then merge all into one file

for (i in seq_along(point_dbs)){
  
  #Work with the sqlite data ----
  
  #connect to the database
  con <- dbConnect(drv=RSQLite::SQLite(), dbname=point_dbs[i])
  print(con)
  
  #list all tables, besides the unnecesary "sqlite_sequence"
  tables <- dbListTables(con)[ dbListTables(con) !="sqlite_sequence"]
  print(tables)
  
  #Initiate the sql_dataframes
  sql_dataframes <- vector("list", length=length(tables))
  
  ## create a dataframe for each table
  for (j in seq_along(tables)) {
    assign(paste0("df",j), dbGetQuery(conn=con, statement=paste("SELECT * FROM '", tables[[j]], "'", sep="")))
  }
  
  #create a shapefile from the dataframes ----
  
  #change the x, y geometry fields to lat/long
  colnames(df2)[colnames(df2) == "ignition_x"]<- "lon"
  colnames(df2)[colnames(df2) == "ignition_y"]<- "lat"
  
  pts_vector <- vect(x = df2,
                     crs = df3$srtext)
  
  #save the file
  writeVector(pts_vector,
           paste0("./", "ignitions_", foa_run),
           filetype= "ESRI Shapefile")
}

