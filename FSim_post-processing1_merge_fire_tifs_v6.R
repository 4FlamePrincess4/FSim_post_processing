library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(RSQLite)

#Set the working directory to the specific outputs folder for the run
wd <- setwd("D:/WFSETP/FSim_Outputs/2022_baseline_t0/okwen_foa2c_r10/")
wd <- setwd("D:/WFSETP/FSim_Outputs/2022_baseline_t0/okwen_foa2c_r10/")

#Create directories to store the results
dir.create("./SeasonFires_merged_tifs/")

#######################################################################################
# NOTE: To run this code, you need to make sure the following FSim outputs are in the #
#       working directory: the FireSizeList.csv files, the  ArrivalDays tifs, and     #
#        the FlameLengths tifs files. You also need to use the rename_tifs script to  #
#       append the directory name (with the run & part info) to the beginning of the  #
#       file names.                                                                   #
#######################################################################################
##                                CODE DESCRIPTION                                   ##
# The plan is to use the fire perimeters to determine whether any fires overlap.      #
# If they do not overlap, we will simply merge the tifs. If they do overlap, we will  #
# store the fire IDs, use the IDs to fetch the correct tifs, merge the perimeters and #
# save the merged extent, clip the tifs to the extent to save processing time, then   #
# do the cell-by-cell minimum arrival day operation to save only the earliest arrival #
# days. We can then resize the tifs back to the foa extent and add them back to the   #
# tif stack in place of the previous version. Lastly, we'll use the perimeters and    #
# ignition points to delete remaining overburn before moving on to the next season.   #
####################################################################################### 

#STEP 1: Record run information below 
###############################################
foa_run <- "FOA2c_r10"
scenario <- ""
run_timepoint <- "baseline_time0"
foa_lcp <- rast("../../../FSim_Run_Files/okwen_foa2c_r1/_inputs/lcp/FOA2c_LCG_LF2022_FBFM40_230_120m.tif")
okwen_perimeter <- st_read("../../../Data/OkWen_shapefiles/FOA_shapefiles/OkWen_AllFOAs_60km_buffer/OkWen_cFOAs_Albers_60km_Buffer.shp")
#foa_extent <- ext(foa_lcp)
#okwen_extent <- ext(okwen_perimeter)
number_of_seasons <- 20000
#Use the below if you have an equal number of seasons for each part
number_of_parts <- 8
seasons_per_part <- c(rep(2500, number_of_parts))
#Use the below alternative if you have different numbers of seasons for each part
seasons_per_part <- c(5000, 7000, 2000, 1000, 5000)

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

#STEP 3: Define a set of functions to use for post-processing.
###############################################

#Create a function to read fire tif files, align them with the foa raster, and stack them
align_and_stack_tifs <- function(file) {
  fire_raster <- rast(file)
  #print(fire_raster)
  aligned_fire <- align_raster(fire_raster, foa_lcp)
  #print(aligned_fire)
  return(aligned_fire)
}
#Create a function to align the individual fire tifs with the FOA lcp tif
align_raster <- function(fire_raster, foa_lcp) {
  # Define the projection of the fire raster to match the projection of the foa raster
  crs(fire_raster) <- crs(foa_lcp)
  # Project the fire raster to align with the foa raster
  fire_raster <- extend(fire_raster, ext(foa_lcp), snap="near")
  return(fire_raster)
}
#Create a function to set non-minimum cell values in a raster stack to NA
set_nonmin_to_na <- function(r_stack) {
  # Replace non-minimum values with NA
  modified_stack <- app(r_stack, fun = function(x) {
    min_vals <- min(x, na.rm = TRUE)
    ifelse(x == min_vals, x, NA)
  })
  
  return(modified_stack)
}
#Create a function to count non-NA values that can be applied to each pixel in a raster stack
count_non_na <- function(x){
  num_non_na <- sum(!is.na(x))
  return(num_non_na)
}
#Create a function to identify which polygon of a multi-polygon layer contains a corresponding ignition point
identify_containing_polygon <- function(polygons, points) {
  sapply(st_geometry(points), function(point) {
    which(sapply(st_geometry(polygons), function(poly) {
      st_contains(st_sfc(poly), st_sfc(point), sparse = FALSE)
    }))
  })
}
#Create a function to delete a polygon containing an ignition point.
delete_polygon_with_point <- function(polygons, point_index) {
  polygons[-point_index, ]
}
#Create a function to use polygons to mask a tif file
mask_tif_with_polygon <- function(raster_stack,polygon, index){
  #Convert the remaining polygon to a SpatVector for terra
  polygon_terra <- vect(polygon)
  #Extract the correct layer from the raster stack using the index
  raster_layer <- raster_stack[[index]]
  #Mask the raster using the polygon
  masked_raster <- mask(raster_layer, polygon_terra, maskvalue=NA)
  return(masked_raster)
}
#Create a function to handle NA values during multiplication
multiply_with_na <- function(x, y) {
  ifel(is.na(x), NA, x * ifel(is.na(y), 0, y))
}
#Create a function to merge a multi-layer spat raster and to keep the first non-NA value in the list, or return an NA if there are no non-NA values.
merge_non_na_first <- function(...) {
  vals <- c(...)
  non_na_vals <- vals[!is.na(vals)]
  if (length(non_na_vals) > 0) {
    return(non_na_vals[1])
  } else {
    return(NA)
  }
}

#STEP 4: Start a for loop to break this up by season and create ArrivalDay and FlameLength
# raster stacks.
########################################################################################################################
# NOTE: Before you can run this code, you need to give the fire tif files unique names indicating their run and part.  #
#   To do this, run the PowerShell script "rename_tifs.ps1", which you can download from the FSim_scripts Google Drive #
#   folder. You need to place the script in the folder, then you can right click and select the option to run the      #
#   script in PowerShell. The script will be renamed as well, but this doesn't matter. I suggest copying the original  #
#   to each folder and running the script before moving on to pasting into the next folder. The script is a small file.#
########################################################################################################################

#Initialize a list to keep track of which fires should not have occurred
# (i.e., ignition inside previous fire perimeter).
fires_to_delete <- list()

for(each_season in unique(firelists$Season)){
  print(paste0("Processing Season ", each_season,"..."))
  this_season_fires <- firelists %>%
    filter(Season == each_season)
  this_season_fireIDs <- this_season_fires$FireID
  this_season_fireIDs <- as.character(this_season_fireIDs)
  this_season_pt <- this_season_fires$Part
  this_season_pt <- as.character(this_season_pt)
  this_season_scen <- this_season_fires$Scenario
  this_season_scen <- as.character(this_season_scen)
  this_season_foa_run <- rep(foa_run, length(this_season_fireIDs))
  #####TODO: Before you run this for real, de-comment the below code and make sure your filenames include the scenario.
  # this_season_AD_filenames <- paste0(wd,"/",this_season_foa_run,"_",this_season_pt,"_ArrivalDays/",
  #                                       this_season_foa_run, "_", this_season_pt, "_",
  #                                       this_season_scen, "_ArrivalDays_FireID_",
  #                                       this_season_fireIDs, ".tif")
  this_season_AD_filenames <- paste0(wd,"/",this_season_foa_run,"_",this_season_pt,"_ArrivalDays/",
                                     this_season_foa_run, "_", this_season_pt, "_ArrivalDays_FireID_",
                                     this_season_fireIDs, ".tif")
  # Read each fire raster, align it with the foa raster, and add it to a raster stack using the align_and_stack_tifs function
  print(paste0("Reading in Arrival Day tifs for Season ", each_season,"..."))
  this_season_AD_stack <- rast(lapply(this_season_AD_filenames, align_and_stack_tifs))
  #Plot the raster stack
  #plot(this_season_AD_stack)
  #print(this_season_AD_stack)
  
  #Now let's list the flame length tif files 
  #####TODO: Before you run this for real, de-comment the below code and make sure your filenames include the scenario.
  # this_season_FL_filenames <- paste0(wd,"/",this_season_foa_run,"_",this_season_pt,"_FlameLengths/",
  #                                    this_season_foa_run, "_", this_season_pt, "_",
  #                                    this_season_scen, "_FlameLengths_FireID_",
  #                                    this_season_fireIDs, ".tif")
  this_season_FL_filenames <- paste0(wd,"/",this_season_foa_run,"_",this_season_pt,"_FlameLengths/",
                                     this_season_foa_run, "_", this_season_pt, "_FlameLengths_FireID_",
                                     this_season_fireIDs, ".tif")
  # Read each fire raster, align it with the foa raster, and add it to a raster stack
  print(paste0("Reading in Flame Length tifs for Season ", each_season,"..."))
  this_season_FL_stack <- rast(lapply(this_season_FL_filenames, align_and_stack_tifs))
  #Plot the raster stack
  #plot(this_season_FL_stack)
  #print(this_season_FL_stack)
  
  #STEP 5: Load the fire perimeters for the season and check whether any overlap.
  ###############################################
  #Load the perimeter database corresponding to the season's part
  con <- dbConnect(SQLite(), dbname = paste0(wd,"/",foa_run, "_",
                                             this_season_pt[1], "_Perimeters.sqlite"))
  # Construct the SQL query to select the polygons based on the IDs
  query1 <- paste("SELECT * FROM perimeters WHERE fire_id IN (", toString(this_season_fireIDs),")")
  # Query the sqlite database to fetch only the perimeters that are in each id list
  season_fire_perims <- dbGetQuery(con, query1)
  #Get the reference system (which is the same for all three)
  ref_sys <- dbGetQuery(con, "SELECT * FROM spatial_ref_sys")
  # Close the database connection
  dbDisconnect(con)
  # Convert the geometry column to sf objects
  season_fire_perims$GEOMETRY <- sf::st_as_sfc(structure(as.list(season_fire_perims$GEOMETRY),class="blob"),
                                             crs = ref_sys$srtext)
  season_fire_perims <- st_as_sf(season_fire_perims)
  # Plot the perimeters
  plot(season_fire_perims, col = "black", max.plot = 1)
  
  #Check whether any of the perimeters overlap
  overlap_matrix <- st_intersects(season_fire_perims, sparse = FALSE)
  #Identify perimeters that overlap with any other perimeters
  overlap_indices <- list()
  n <- nrow(overlap_matrix)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (overlap_matrix[i, j] && i != j) {
        overlap_indices <- append(overlap_indices, list(c(i, j)))
      }
    }
  }
 
  
  #Print how many fires overlap in the season.
  print(paste0("There are ",length(overlap_indices), " cases where fires overlap in season ", each_season, "."))
  
  #STEP 6: If there is no overlap, merge the fire tifs and write to the folders.
  ###############################################
  #If no fires overlap
  if(length(overlap_indices) == 0){
    #First, create a raster that identifies each fire in each pixel
    # Create a copy of the raster stack to store the fire ID values
    this_season_ID_stack <- rast(nrows = nrow(this_season_AD_stack), 
                                    ncols = ncol(this_season_AD_stack), 
                                    ext = ext(this_season_AD_stack), 
                                    crs = crs(this_season_AD_stack), 
                                    vals = NA)
    this_season_fireIDs <- as.numeric(this_season_fireIDs)
    print("Creating the fire ID raster...")
    for(fireid in 1:nlyr(this_season_AD_stack)){
      this_fireid <- this_season_fireIDs[fireid]
      this_AD_raster <- this_season_AD_stack[[fireid]]
      # Replace non-NA values with the fire ID directly in the final raster 
      #  with a mask to avoid memory issues.
      mask <- !is.na(this_AD_raster)
      this_season_ID_stack[mask] <- this_fireid
    }
    print("Merging the AD and FL rasters...")
    #Merge the ArrivalDay raster stack
    season_fires_AD_raster <- app(this_season_AD_stack, fun = merge_non_na_first)
    #Merge the FlameLength raster stack
    season_fires_FL_raster <- app(this_season_FL_stack, fun = merge_non_na_first)
    
    #Create a stack of the three rasters
    season_fires_raster_stack <- c(this_season_ID_stack, season_fires_AD_raster, 
                                   season_fires_FL_raster)
    names(season_fires_raster_stack) <- c("Fire_IDs", "Julian_Arrival_Days", "Flame_Lengths_ft")
    
    #Plot to check the stack merged correctly
    plot(season_fires_raster_stack, main = paste0("Season ", each_season))
    
    
    #Write the merged season arrival day raster
    writeRaster(season_fires_raster_stack, filename=paste0("./SeasonFires_merged_tifs/Season",
                                                        each_season,"_merged_IDs_ADs_FLs.tif"),
                overwrite = TRUE)

  }else{ #End of if-else scenario where the fire perimeters do not overlap.
    
    #STEP 7: If there is overlap, use the extent of the merged perimeters to clip the tifs,
    # then keep only the minimum arrival day values. Use these to filter the flame length tifs.
    ###############################################
    #First, we need to identify the fires that overlap
    # Iterate over overlap_indices to extract fire IDs
    overlapping_fire_ids <- list()
    for (pair in overlap_indices) {
      fire_id1 <- season_fire_perims$fire_id[pair[1]]
      fire_id2 <- season_fire_perims$fire_id[pair[2]]
      # Store fire IDs in the list
      overlapping_fire_ids <- append(overlapping_fire_ids, list(data.frame(fire_id1 = fire_id1, fire_id2 = fire_id2)))
    }
    # Combine the list of data frames into a single data frame
    overlapping_fire_ids_df <- do.call(rbind, overlapping_fire_ids)
    overlapping_fire_indices_df <- do.call(rbind, overlap_indices)
    overlapping_fire_indices_df <- as.data.frame(overlapping_fire_indices_df)
    colnames(overlapping_fire_indices_df) <- c("fire_index1", "fire_index2")
    
    #Next, we need to check whether there are any cases where more than two fires overlap 
    # at a location on a pixel-by-pixel basis
    # Compute the number of non-NA values per pixel
    print("Checking the number of non-NA values (# of fires) per pixel...")
    num_non_na_per_pixel <- app(this_season_AD_stack, fun=count_non_na)
    max_overlapping_fires <- global(num_non_na_per_pixel, max, na.rm=TRUE)[[1]]
    if(max_overlapping_fires > 2){
      print("There are more than two overlapping fires at a single pixel.")
      # More than two overlapping fires at this pixel
      cells_over_two_overlap <- which(values(num_non_na_per_pixel) > 2)
      excess_overlap_fires <- c()
      for(i in seq_along(cells_over_two_overlap)){
        # Extract arrival days for all fires at the pixel where the excess overlap occurs
        arrival_days_df <- extract(this_season_AD_stack, i)
        # Convert the data frame to a vector
        arrival_days_vector <- as.vector(t(arrival_days_df[1,]))
        # Find the indices of non-NA values
        non_na_indices <- which(!is.na(arrival_days_vector))
        excess_overlap_fireids <- overlapping_fire_ids[non_na_indices]
        excess_overlap_fires <- c(excess_overlap_fires, excess_overlap_fireids)
      }
      excess_fires <- unique(excess_overlap_fires)
      print(paste0("These are the excess fires: ", excess_fires, ". Consider revising the code to deal with this case."))
      
      #### TODO: Handle the case with more than two overlapping fires (see v2 of this code). ####
      
    }else if(max_overlapping_fires > 1 && max_overlapping_fires <=2){ #End of if-else scenario where there are more than two overlapping fires somewhere on the landscape.
      print("There are no more than two fires overlapping at any given pixel on the landscape.")
      # There are only two overlapping fires in each case of overlap.
      #Print the overlapping fire IDs
      print(paste0("These are the cases: "))
      print(overlapping_fire_ids_df)
      
      #Read in ignitions for the overlapping fire ids
      #We only need the unique ids for this
      unique_overlapping_fire_ids <- c(overlapping_fire_ids_df$fire_id1, 
                                       overlapping_fire_ids_df$fire_id2)
      unique_overlapping_fire_ids <- unique(unique_overlapping_fire_ids)
      # Load the ignition database corresponding to the season part
      con <- dbConnect(SQLite(), dbname = paste0(wd,"/",foa_run, "_",
                                                 this_season_pt[1], "_Ignitions.sqlite"))
      # Construct the SQL query to select the ignitions based on the IDs
      query <- paste("SELECT * FROM ignitions WHERE fire_id IN (", 
                     toString(unique_overlapping_fire_ids),")")
      # Query the sqlite database to fetch only the ignitions that are in the overlapping id list
      overlapping_ig <- dbGetQuery(con, query)
      #Get the reference system (which is the same for all three)
      ref_sys <- dbGetQuery(con, "SELECT * FROM spatial_ref_sys")
      # Close the database connection
      dbDisconnect(con)
      # Convert the geometry column to sf objects
      colnames(overlapping_ig)[colnames(overlapping_ig) == "ignition_x"]<- "lon"
      colnames(overlapping_ig)[colnames(overlapping_ig) == "ignition_y"]<- "lat"
      overlapping_ig <- terra::vect(x = overlapping_ig, crs=ref_sys$srtext)
      overlapping_ig <- as_sf(overlapping_ig)
      
      #Use the indices to subset the ArrivalDay and FlameLength stacks
      #We don't need the pairs of fires for this operation. The arrival day info will guide which values to keep.
      unique_overlapping_fire_indices <- c(overlapping_fire_indices_df$fire_index1, 
                                           overlapping_fire_indices_df$fire_index2)
      unique_overlapping_fire_indices <- unique(unique_overlapping_fire_indices)
      overlapping_AD_stack <- this_season_AD_stack[[unique_overlapping_fire_indices]]
      overlapping_FL_stack <- this_season_FL_stack[[unique_overlapping_fire_indices]]
      #Now perform the minimum arrival day operation on the AD stack
      #Find the minimum of each cell across a stack of fire rasters and set all other values to NA.
      print(paste0("Setting non-minimum arrival days of overlapping fires to NA."))
      #earliest_arrival_AD_stack <- set_nonmin_to_na(overlapping_AD_stack)
      #Let's try a different approach to see if it's faster
      # Create a mask for pixels with exactly 2 non-NA values
      non_na_mask <- (num_non_na_per_pixel == 2)
            # Identify minimum values in the stack
      min_vals_stack <- app(this_season_AD_stack, fun = function(x) {
        min_val <- min(x, na.rm = TRUE)
        return(ifelse(x == min_val, x, NA))
      })
            # Apply the non_na_mask to keep only the relevant pixels
      earliest_arrival_AD_stack <- mask(min_vals_stack, non_na_mask, maskvalues = FALSE)
      
      #Mask the FL rasters with the edited AD rasters to keep only FL values for the earliest arrival days.
      print(paste0("Masking the flame length rasters to set non-min arrival days of overlapping fires to NA."))
      overlapping_FL_stack <- align_raster(overlapping_FL_stack, earliest_arrival_AD_stack)
      earliest_arrival_FL_stack <- mask(overlapping_FL_stack, earliest_arrival_AD_stack)
      plot(earliest_arrival_FL_stack, 
           main = paste0("Corrected flame length rasters for overlapping fires: ", overlapping_fire_ids))
      
      
      #STEP 8: Now that we've deleted overburn on a pixel-by-pixel basis using the arrival day information,
      # we still need to delete any fires that should not have ignited within the burn scar of previous fires
      # We will also plot any cases of overlapping fires in order to devise code to handle any weird cases that arise.
      ###############################################
      
      #First, check for cases where the ignitions fall into the corresponding polygons.
      # Iterate over each pair of overlapping fires.
      for(pairs in 1:nrow(overlapping_fire_ids_df)){
        # Grab each ID
        first_ID <- overlapping_fire_ids_df$fire_id1[pairs]
        second_ID <- overlapping_fire_ids_df$fire_id2[pairs]
        # Grab each index
        first_index <- overlapping_fire_indices_df$fire_index1[pairs]
        second_index <- overlapping_fire_indices_df$fire_index2[pairs]
        # Subset the season polygons and the overlapping fire ignitions with these IDs
        first_ig <- overlapping_ig %>%
          dplyr::filter(fire_id == first_ID)
        second_ig <- overlapping_ig %>%
          dplyr::filter(fire_id == second_ID)
        first_perim <- season_fire_perims %>%
          dplyr::filter(fire_id == first_ID)
        second_perim <- season_fire_perims %>%
          dplyr::filter(fire_id == second_ID)
        # Plot these as a sanity check
        overlap_case_plot <- ggplot()+
          geom_sf(data=first_perim, aes(color = "Fire 1 Perim"))+
          geom_sf(data=first_ig, aes(color = "Fire 1 Ig"))+
          geom_sf(data=second_perim, aes(color = "Fire 2 Perim"))+
          geom_sf(data=second_ig, aes(color = "Fire 2 Ig"))+
          scale_color_manual(values=c("Fire 1 Perim" = "black","Fire 1 Ig" = "red",
                                      "Fire 2 Perim" = "blue","Fire 2 Ig" = "green"))+
          theme_minimal()
        print(overlap_case_plot)
        
    # Check whether the first ignition falls in the second perimeter
    first_ignition_in_second_perim <- st_intersects(first_ig, second_perim, sparse = FALSE)
    # Check whether the second ignition falls in the first perimeter
    second_ignition_in_first_perim <- st_intersects(second_ig, first_perim, sparse = FALSE)
    #If the first ignition is inside the second perimeter
    if(first_ignition_in_second_perim){
      #Check whether the ignition day is earlier than the arrival day at the pixel of the perimeter
      # Extract the coordinates of the ignition point
      first_ig_coords <- st_coordinates(first_ig)
      # Extract the arrival day value at the ignition point location
      second_fire_AD <- this_season_AD_stack[[second_index]]
      second_fire_AD_df <- extract(second_fire_AD, matrix(first_ig_coords, ncol=2))
      # Extract the value itself from the resulting dataframe
      second_fire_AD_value <- second_fire_AD_df[1,1]
      # Compare the start_day with the arrival_day_value
      if (!is.na(second_fire_AD_value)) {
        if (first_fire_ig$start_day <= second_fire_AD_value) {
          print("The ignition date for fire 1 is less than or equal to the arrival day value for fire 2 at the ignition point's location.")
          print("We will keep both fires.")
          } else { #End of if-else scenario where the first fire ignition is within the second fire perim
              # & the first fire would have prevented the second fire from spreading. #End of if-else scenario where the ignition date of fire 1 is less than or equal to the AD of fire 2.
            print("The ignition date for fire 1 is greater than the arrival day value for fire 2 at the ignition point's location.")
            print("Fire 1 should not have occurred and will be deleted.")
            #Store the ID in a list to delete the FL and AD rasters from the stacks just before merging.
            # Wait to do this to avoid messing up the stack indexing.
            fires_to_delete <- append(fires_to_delete, first_index)
          } #End of if-else scenario where the ignition for fire 1 is inside of the fire 2 perimeter and
          # fire 1 should not have ignited.
          } else { #End of scenario where the second fire AD value is not NA.
            print("No valid arrival_day value at the ignition point's location.")
            print(paste0("Ignition fire ID = ", first_ID, ". Perimeter fire ID = ", second_ID, "."))
          } #End of scenario where the second fire AD value is NA.
      } else if(second_ignition_in_first_perim){ #End of if-else scenario where the first fire ignition is inside of the second fire perimeter.
        #Check whether the ignition day is earlier than the arrival day at the pixel of the perimeter
        # Extract the coordinates of the ignition point
        second_ig_coords <- st_coordinates(second_ig)
        # Extract the arrival day value at the ignition point location
        first_fire_AD <- this_season_AD_stack[[first_index]]
        first_fire_AD_df <- extract(first_fire_AD, matrix(second_ig_coords, ncol=2))
        # Extract the value itself from the resulting dataframe
        first_fire_AD_value <- first_fire_AD_df[1,1]
        # Compare the start_day with the arrival_day_value
        if (!is.na(first_fire_AD_value)) {
          if (second_ig$start_day <= first_fire_AD_value) {
            print("The ignition date for fire 2 is less than or equal to the arrival day value for fire 1 at the ignition point's location.")
            print("We will keep both fires.")
            } else { #End of if-else scenario where the fire 2 ignition is within the fire 1 perim
              # & fire 2 would have prevented fire 1 from spreading.
            #End of if-else scenario where the fire 2 ignition is earlier than the fire 1 AD value.
            print("The ignition date for fire 2 is greater than the arrival day value for fire 1 at the ignition point's location.")
            print("Fire 2 should not have occurred and will be deleted.")
            #Store the ID in a list to delete the FL and AD rasters from the stacks just before merging.
            # Wait to do this to avoid messing up the stack indexing.
            fires_to_delete <- append(fires_to_delete, second_index)
            } #End of if-else scenario where the ignition for fire 1 is inside of the fire 2 perimeter and the
            # ignition date is later than the fire 2 arrival day, so fire 1 would not have occurred.
          } else {#End of scenario where the AD value at the ignition point is not NA.
            print("No valid arrival_day value at the ignition point's location.")
            print(paste0("Ignition fire ID = ", second_ID, ". Perimeter fire ID = ", first_ID, "."))
            } #End of if-else scenario where the AD value at the ignition point is NA.
        } else if(!(first_ignition_in_second_perim) && !(second_ignition_in_first_perim)){ #End of scenario where fire 2 ignition is inside of fire 1 perimeter.
          #If the ignitions don't fall into the other polygons,
          # subtract the polygon of the ignition point of fire 1 from the polygon of the later fire (fire 2)
          print("The ignitions don't fall into the perimeters of the overlapping fires.")
          print("We are assuming that keeping only the earliest arrival day value will take care of overburn.")
          } #End of the if-else situation where the ignitions don't overlap with the polygons of overlapping fires.
        } #End of the for loop to edit each pair of fires that overlap.
      } else { #End of if-else scenario where there are only two overlapping fires in each case of overlap.
        print("The polygons overlapped but the raster values do not overlap. We will set non-min arrival day values to NA to be sure.")
        #Use the indices to subset the ArrivalDay and FlameLength stacks
        #We don't need the pairs of fires for this operation. The arrival day info will guide which values to keep.
        unique_overlapping_fire_indices <- c(overlapping_fire_indices_df$fire_index1,
                                             overlapping_fire_indices_df$fire_index2)
        unique_overlapping_fire_indices <- unique(unique_overlapping_fire_indices)
        overlapping_AD_stack <- this_season_AD_stack[[unique_overlapping_fire_indices]]
        overlapping_FL_stack <- this_season_AD_stack[[unique_overlapping_fire_indices]]
        #Mask any non-minimum values with NA
        # Create a mask for pixels with exactly 2 non-NA values
        non_na_mask <- (num_non_na_per_pixel == 2)
        # Identify minimum values in the stack
        min_vals_stack <- app(this_season_AD_stack, fun = function(x) {
          min_val <- min(x, na.rm = TRUE)
          return(ifelse(x == min_val, x, NA))
        })
        # Apply the non_na_mask to keep only the relevant pixels
        earliest_arrival_AD_stack <- mask(min_vals_stack, non_na_mask, maskvalues = FALSE)
        #Mask the FL rasters with the edited AD rasters to keep only FL values for the earliest arrival days.
        print(paste0("Masking the flame length rasters to set non-min arrival days of overlapping fires to NA."))
        overlapping_FL_stack <- align_raster(overlapping_FL_stack, earliest_arrival_AD_stack)
        earliest_arrival_AD_stack <- mask(overlapping_FL_stack, earliest_arrival_AD_stack)
        plot(earliest_arrival_FL_stack,
             main = paste0("Corrected flame length rasters for overlapping fires: ", overlapping_fire_ids))
      }
    #Overwrite the corresponding tifs in the season AD & FL stacks
    for(i in seq_along(unique_overlapping_fire_indices)){
      this_AD <- earliest_arrival_AD_stack[[i]]
      this_FL <- earliest_arrival_FL_stack[[i]]
      this_index <- unique_overlapping_fire_indices[i]
      print(paste0("Replacing fire ", this_index, " in the original stack with the revised version."))
      this_season_AD_stack[[this_index]] <- this_AD
      rm(this_AD)
      this_season_FL_stack[[this_index]] <- this_FL
      # Remove temporary objects and run garbage collection
      rm(this_FL)
      gc()
    } #End of for loop to replace AD & FL rasters with edited AD & FL rasters
    
    #Delete the arrival day and flame length tifs of fires that shouldn't have happened
     if(length(fires_to_delete) > 0){
       fires_to_delete <- unlist(fires_to_delete)
       earliest_arrival_fire_stack <- earliest_arrival_fire_stack[-fires_to_delete]
       earliest_arrival_FL_stack <- earliest_arrival_FL_stack[-fires_to_delete]
     }

    #Create a raster that identifies each fire in each pixel
    # Create a copy of the raster stack to store the fire ID values
    this_season_ID_stack <- rast(nrows = nrow(this_season_AD_stack), 
                                 ncols = ncol(this_season_AD_stack), 
                                 ext = ext(this_season_AD_stack), 
                                 crs = crs(this_season_AD_stack), 
                                 vals = NA)
    this_season_fireIDs <- as.numeric(this_season_fireIDs)
    print("Creating the fire ID raster...")
    for(fireid in 1:nlyr(this_season_AD_stack)){
      this_fireid <- this_season_fireIDs[fireid]
      this_AD_raster <- this_season_AD_stack[[fireid]]
      # Replace non-NA values with the fire ID directly in the final raster 
      #  with a mask to avoid memory issues.
      mask <- !is.na(this_AD_raster)
      this_season_ID_stack[mask] <- this_fireid
    }
    print("Merging the AD and FL rasters...")
    #Merge the ArrivalDay raster stack
    season_fires_AD_raster <- app(this_season_AD_stack, fun = merge_non_na_first)
    #Merge the FlameLength raster stack
    season_fires_FL_raster <- app(this_season_FL_stack, fun = merge_non_na_first)
    
    #Create a stack of the three rasters
    season_fires_raster_stack <- c(this_season_ID_stack, season_fires_AD_raster, 
                                   season_fires_FL_raster)
    names(season_fires_raster_stack) <- c("Fire_IDs", "Julian_Arrival_Days", "Flame_Lengths_ft")
    
    #Plot to check the stack merged correctly
    plot(season_fires_raster_stack, main = paste0("Season ", each_season))
    
    
    #Write the merged season arrival day raster
    writeRaster(season_fires_raster_stack, filename=paste0("./SeasonFires_merged_tifs/Season",
                                                           each_season,"_merged_IDs_ADs_FLs.tif"),
                overwrite = TRUE)
  } #End of the if-else scenario where fires overlap.
} #End of for loop to merge AD & FL tifs for each season.
