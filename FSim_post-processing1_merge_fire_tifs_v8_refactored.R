library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(RSQLite)
library(future)
library(furrr)

#Set the working directory to the specific outputs folder for the run
wd <- setwd("C:/FSimProjectFolder/FSim_Outputs/okwen_2022_baseline_foa3c_r3_time0/")
wd <- setwd("C:/FSimProjectFolder/FSim_Outputs/okwen_2022_baseline_foa3c_r3_time0/")

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

#STEP 3: Define a set of helper functions 
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
set_nonmin_to_na <- function(r_stack, non_na_mask) {
  # Replace non-minimum values with NA only where non_na_mask is TRUE
  modified_stack <- app(r_stack, fun = function(x, mask) {
    if (mask) {
      min_vals <- min(x, na.rm = TRUE)
      return(ifelse(x == min_vals, x, NA))
    } else {
      return(x)
    }
  }, non_na_mask)
  
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
merge_non_na_first <- function(...){
  vals <- c(...)
  non_na_vals <- vals[!is.na(vals)]
  if(length(non_na_vals) > 0){
    return(non_na_vals[1])
  }else{
    return(NA)
  }
}
#Create a function to find overlap indices
find_overlap_indices <- function(overlap_matrix) {
  overlap_indices <- list()
  n <- nrow(overlap_matrix)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (overlap_matrix[i, j] && i != j) {
        overlap_indices <- append(overlap_indices, list(c(i, j)))
      }
    }
  }
  return(overlap_indices)
}

#STEP 4: Define smaller functions to be used within the larger process_fire_season function, which is then run with future_map()
########################################################################################################################
# NOTE: Before you can run this code, you need to give the fire tif files unique names indicating their run and part.  #
#   To do this, run the PowerShell script "rename_tifs.ps1", which you can download from the FSim_scripts Google Drive #
#   folder. You need to place the script in the folder, then you can right click and select the option to run the      #
#   script in PowerShell. The script will be renamed as well, but this doesn't matter. I suggest copying the original  #
#   to each folder and running the script before moving on to pasting into the next folder. The script is a small file.#
########################################################################################################################

#Define a function to process a single season
process_single_fire_season <- function(this_season_AD_stack, this_season_FL_stack, each_season, wd, this_season_fireIDs) {
  #This function processes a season where there is only a single fire.
  print(paste0("There is only one fire in season ", each_season))
  #Create the fire ID raster stack
  this_season_ID_stack <- rast(nrows = nrow(this_season_AD_stack), ncols = ncol(this_season_AD_stack), ext = ext(this_season_AD_stack), crs = crs(this_season_AD_stack), vals = NA)
  this_season_fireIDs <- as.numeric(this_season_fireIDs)
  print("Creating the fire ID raster...")
  for(fireid in 1:nlyr(this_season_AD_stack)){
    this_fireid <- this_season_fireIDs[fireid]
    this_AD_raster <- this_season_AD_stack[[fireid]]
    mask <- !is.na(this_AD_raster)
    this_season_ID_stack[mask] <- this_fireid
  }
  #Combine the fire ID stack with the AD & FL stacks. (Each "stack" actually only has one layer because there is only one fire.)
  season_fires_raster_stack <- c(this_season_ID_stack, this_season_AD_stack, this_season_FL_stack)
  names(season_fires_raster_stack) <- c("Fire_IDs", "Julian_Arrival_Days", "Flame_Lengths_ft")
  plot(season_fires_raster_stack, main = paste0("Season ", each_season))
  #Write the resulting 3-band raster stack.
  writeRaster(season_fires_raster_stack, filename=paste0("./SeasonFires_merged_tifs/Season", each_season,"_merged_IDs_ADs_FLs.tif"), overwrite = TRUE)
}

merge_and_write_rasters <- function(this_season_AD_stack, this_season_FL_stack, each_season, wd, this_season_fireIDs) {
  #First, make a fire ID stack to identify each pixel burned by a given fire.
  this_season_ID_stack <- rast(nrows = nrow(this_season_AD_stack), ncols = ncol(this_season_AD_stack), ext = ext(this_season_AD_stack), crs = crs(this_season_AD_stack), vals = NA)
  this_season_fireIDs <- as.numeric(this_season_fireIDs)
  print("Creating the fire ID raster...")
  for(fireid in 1:nlyr(this_season_AD_stack)){
    this_fireid <- this_season_fireIDs[fireid]
    this_AD_raster <- this_season_AD_stack[[fireid]]
    mask <- !is.na(this_AD_raster)
    this_season_ID_stack[mask] <- this_fireid
  }
  #Next, merge the AD and FL rasters with the helper function merge_non_na_first
  print("Merging the AD and FL rasters...")
  season_fires_AD_raster <- app(this_season_AD_stack, fun = merge_non_na_first)
  season_fires_FL_raster <- app(this_season_FL_stack, fun = merge_non_na_first)
  #Create the 3-band season raster stack
  season_fires_raster_stack <- c(this_season_ID_stack, season_fires_AD_raster, season_fires_FL_raster)
  names(season_fires_raster_stack) <- c("Fire_IDs", "Julian_Arrival_Days", "Flame_Lengths_ft")
  plot(season_fires_raster_stack, main = paste0("Season ", each_season, " fire IDs, ADs, & FLs"))
  #Save the 3-band season raster stack
  writeRaster(season_fires_raster_stack, filename=paste0("./SeasonFires_merged_tifs/Season", each_season,"_merged_IDs_ADs_FLs.tif"), overwrite = TRUE)
}

handle_more_than_two_overlaps <- function(this_season_AD_stack, this_season_FL_stack, each_season, wd, overlap_indices, season_fire_perims, ref_sys, firelists, foa_run) {
  #This function processes cases where there are more than two overlapping fires at a single pixel.
  print(paste("There are up to", max_overlapping_fires, "fires at a single pixel."))
  #Find which cell indices have the excess overlap
  cells_over_two_overlap <- which(values(num_non_na_per_pixel) > 2)
  excess_overlap_fires <- c()
  non_na_indices_list <- list()
  #Get the fire IDs for the fires that overlap at these locations
  for(i in seq_along(cells_over_two_overlap)){
    # Extract arrival days for all fires at the pixel where the excess overlap occurs
    arrival_days_df <- extract(this_season_AD_stack, cells_over_two_overlap[i])
    # Convert the data frame to a vector
    arrival_days_vector <- as.vector(t(arrival_days_df[1,]))
    # Find the indices of non-NA values
    non_na_indices <- which(!is.na(arrival_days_vector))
    excess_overlap_fireids <- this_season_fireIDs[non_na_indices]
    excess_overlap_fires <- c(excess_overlap_fires, excess_overlap_fireids)
    non_na_indices_list[[i]] <- non_na_indices
  }
  #Simplify to the unique IDs for these fires and report out
  excess_fires <- unique(excess_overlap_fires)
  print(paste("These are the excess fires:", paste(excess_fires, collapse = ", ")))
  
  # Load the ignitions and perimeters
  # Load the ignition database corresponding to the season part
  con <- dbConnect(RSQLite::SQLite(), dbname = paste0(wd,"/",foa_run, "_", this_season_pt[1], "_Ignitions.sqlite"))
  # Construct the SQL query to select the ignitions based on the IDs
  query <- paste("SELECT * FROM ignitions WHERE fire_id IN (", toString(excess_fires),")")
  # Query the sqlite database to fetch only the ignitions that are in the overlapping id list
  overlapping_ig <- dbGetQuery(con, query)
  # Get the reference system (which is the same for all three)
  ref_sys <- dbGetQuery(con, "SELECT * FROM spatial_ref_sys")
  # Close the database connection
  dbDisconnect(con)
  # Convert the geometry column to sf objects
  colnames(overlapping_ig)[colnames(overlapping_ig) == "ignition_x"]<- "lon"
  colnames(overlapping_ig)[colnames(overlapping_ig) == "ignition_y"]<- "lat"
  overlapping_ig <- terra::vect(x = overlapping_ig, crs=ref_sys$srtext)
  overlapping_ig <- as_sf(overlapping_ig)
  
  # Use the indices to subset the ArrivalDay and FlameLength stacks
  overlapping_AD_stack <- this_season_AD_stack[[unique(unlist(non_na_indices_list))]]
  overlapping_FL_stack <- this_season_FL_stack[[unique(unlist(non_na_indices_list))]]
  # Find the minimum of each cell across a stack of fire rasters and set all other values to NA.
  print("Setting non-minimum arrival days of overlapping fires to NA.")
  # Create a mask for pixels with 2 or more non-NA values to improve efficiency
  non_na_mask <- num_non_na_per_pixel >= 2
  # Identify minimum values in the stack
  min_vals_stack <- app(this_season_AD_stack, fun = function(x){
    min_val <- min(x, na.rm = TRUE)
    return(ifelse(x == min_val, x, NA))
  })
  #Use these to mask the arrival day stack, setting non-minimum values to NA
  earliest_arrival_AD_stack <- mask(min_vals_stack, non_na_mask, maskvalues = FALSE)
  plot(earliest_arrival_AD_stack, main = paste("Corrected arrival day rasters for overlapping fires:", paste(excess_fires, collapse = ", ")))
  # Mask the FL rasters with the edited AD rasters to keep only FL values for the earliest arrival days.
  print("Masking the flame length rasters to set non-min arrival days of overlapping fires to NA.")
  overlapping_FL_stack <- align_raster(overlapping_FL_stack, earliest_arrival_AD_stack)
  earliest_arrival_FL_stack <- mask(overlapping_FL_stack, earliest_arrival_AD_stack)
  plot(earliest_arrival_FL_stack, main = paste("Corrected flame length rasters for overlapping fires:", paste(excess_fires, collapse = ", ")))
  
  # The non-minimum values have been set to NA for more than three fires, but we still need to deal 
  # with cases where the ignitions of fires
  # Grab just the overlapping perimeters
  overlapping_perims <- season_fire_perims[season_fire_perims$fire_id %in% excess_fires, ]
  fires_to_delete <- data.frame(fire_id = character(), stringsAsFactors = FALSE)
  #Check for ignitions within each overlapping perimeter
  for (i in seq_along(overlapping_perims)) {
    # Get the current perimeter and raster
    perimeter <- overlapping_perims[i, ]
    arrival_raster <- earliest_arrival_AD_stack[[i]]
    # Find ignition points within this perimeter
    points_in_perimeter <- terra::intersect(overlapping_ig, perimeter)
    # If there are no points in the perimeter, skip to the next iteration
    if (nrow(points_in_perimeter) == 0) next
    # Loop through the points within the perimeter
    for (j in 1:nrow(points_in_perimeter)) {
      # Get the ignition point and its start day
      ignition_point <- points_in_perimeter[j, ]
      start_day <- ignition_point$start_day
      ignition_id <- ignition_point$fire_id
      ignition_index <- which(this_season_fireIDs == ignition_id)
      # Extract the arrival day at the ignition point location
      arrival_day <- terra::extract(arrival_raster, terra::geom(ignition_point))[[1]]
      # Compare start day with arrival day
      if (start_day > arrival_day) {
        fires_to_delete <- rbind(fires_to_delete, data.frame(fire_id = ignition_id))
      }
    
    # Plot to see what the cases look like
    ggplot() +
      geom_sf(data = overlapping_perims, aes(fill = as.factor(fire_id)), alpha = 0.5, color = "black") +
      geom_sf_text(data = overlapping_perims, aes(label = fire_id), size = 5, color = "blue") +
      geom_sf(data = overlapping_ig, aes(color = as.factor(fire_id)), size = 3) +
      geom_sf_text(data = overlapping_ig, aes(label = fire_id), nudge_y = 0.05, color = "red") +
      theme_minimal() +
      labs(title = "Fire Perimeters and Ignitions", fill = "Fire ID", color = "Ignition Fire ID")
    }
  }
}


handle_two_or_fewer_overlaps <- function(this_season_AD_stack, this_season_FL_stack, overlapping_fire_ids_df, overlapping_fire_indices_df, num_non_na_per_pixel, season_fire_perims, wd, firelists, foa_run) {
  #If there are only two overlaps
  print("There are at most two fires overlapping at any given pixel.")
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
  # Create a mask for pixels with exactly 2 non-NA values to improve efficiency
  non_na_mask <- num_non_na_per_pixel == 2
  #earliest_arrival_AD_stack <- set_nonmin_to_na(overlapping_AD_stack, non_na_mask)
  #The terra app functiond can't accept a SpatRaster as a second argument, so let's try a different approach.
  # Identify minimum values in the stack
  min_vals_stack <- app(this_season_AD_stack, fun = function(x){
    min_val <- min(x, na.rm = TRUE)
    return(ifelse(x == min_val, x, NA))
  })
  earliest_arrival_AD_stack <- mask(min_vals_stack, non_na_mask, maskvalues = FALSE)
  #Mask the FL rasters with the edited AD rasters to keep only FL values for the earliest arrival days.
  print(paste0("Masking the flame length rasters to set non-min arrival days of overlapping fires to NA."))
  overlapping_FL_stack <- align_raster(overlapping_FL_stack, earliest_arrival_AD_stack)
  earliest_arrival_FL_stack <- mask(overlapping_FL_stack, earliest_arrival_AD_stack)
  plot(earliest_arrival_FL_stack, 
       main = paste0("Corrected flame length rasters for overlapping fires: ", overlapping_fire_ids))
}

process_overlapping_fires <- function(this_season_AD_stack, this_season_FL_stack, each_season, wd, overlap_indices, season_fire_perims, ref_sys, firelists, foa_run, align_raster, mask, app, extract, global, dbConnect, dbGetQuery, dbDisconnect) {
  #Create a dataframe with overlapping fire IDs
  overlapping_fire_ids_df <- do.call(rbind, lapply(overlap_indices, function(pair) {
    data.frame(fire_id1 = season_fire_perims$fire_id[pair[1]], fire_id2 = season_fire_perims$fire_id[pair[2]])
  }))
  #Create a corresponding dataframe of overlapping fire indices
  overlapping_fire_indices_df <- do.call(rbind, overlap_indices)
  colnames(overlapping_fire_indices_df) <- c("fire_index1", "fire_index2")
  #Check how many non-na values there are per pixel to determine which overlap function to apply.
  print("Checking the number of non-NA values (# of fires) per pixel...")
  num_non_na_per_pixel <- app(this_season_AD_stack, fun=count_non_na)
  max_overlapping_fires <- global(num_non_na_per_pixel, max, na.rm=TRUE)[[1]]
  #If there are more than two overlapping fires, use the handle_more_than_two_overlaps function.
  if(max_overlapping_fires > 2) {
    handle_more_than_two_overlaps(this_season_AD_stack, this_season_FL_stack, each_season, wd, overlap_indices, season_fire_perims, ref_sys, firelists, foa_run)
  } else { #If there are only two fires at any given pixel, use the handle_two_or_fewer_overlaps function.
    handle_two_or_fewer_overlaps(this_season_AD_stack, this_season_FL_stack, overlapping_fire_ids_df, overlapping_fire_indices_df, num_non_na_per_pixel, season_fire_perims, wd, firelists, foa_run)
  }
  #Overwrite the corresponding tifs in the season AD & FL stacks
  for(i in seq_along(unique_overlapping_fire_indices)){
    this_AD <- earliest_arrival_AD_stack[[i]]
    this_FL <- earliest_arrival_FL_stack[[i]]
    this_index <- unique_overlapping_fire_indices[i]
    this_ID <- unique_overlapping_fire_ids[i]
    print(paste0("Replacing fire ", this_ID, " in the original stack with the revised version."))
    this_season_AD_stack[[this_index]] <- this_AD
    rm(this_AD)
    this_season_FL_stack[[this_index]] <- this_FL
    # Remove temporary objects and run garbage collection
    rm(this_FL)
    gc()
  }
  #Delete the arrival day and flame length tifs of fires that shouldn't have happened
  if(length(fires_to_delete) > 0){
    fires_to_delete <- unlist(fires_to_delete)
    this_season_AD_stack <- this_season_AD_stack[[-fires_to_delete]]
    this_season_FL_stack <- this_season_FL_stack[[-fires_to_delete]]
    this_season_fireIDs <- this_season_fireIDs[-fires_to_delete]
  }
  #Finally, export the resulting edited tifs in a 3-band season raster stack
  merge_and_write_rasters(this_season_AD_stack, this_season_FL_stack, each_season, wd, this_season_fireIDs)
}

process_fire_season <- function(each_season, firelists, wd, foa_run, align_and_stack_tifs, dbConnect, dbGetQuery, dbDisconnect, merge_non_na_first, count_non_na, app, global, mask, extract, align_raster, writeRaster) {
  print(paste0("Processing Season ", each_season,"..."))
  #Subset the firelists by the current season
  this_season_fires <- firelists %>%
    filter(Season == each_season)
  #Fetch vectors of other run information
  this_season_fireIDs <- as.character(this_season_fires$FireID)
  this_season_pt <- as.character(this_season_fires$Part)
  this_season_scen <- as.character(this_season_fires$Scenario)
  this_season_foa_run <- rep(foa_run, length(this_season_fireIDs))
  #Use this info to read in the arrival day tif filenames for this season's fires
  this_season_AD_filenames <- paste0(wd,"/",this_season_foa_run,"_",this_season_pt,"_ArrivalDays/",
                                     this_season_foa_run, "_", this_season_pt, "_", this_season_scen, "_ArrivalDays_FireID_",
                                     this_season_fireIDs, ".tif")
  print(paste0("Reading in Arrival Day tifs for Season ", each_season,"..."))
  this_season_AD_stack <- rast(lapply(this_season_AD_filenames, align_and_stack_tifs))
  #Use this info to read in the flame length tif filenames for this season's fires
  this_season_FL_filenames <- paste0(wd,"/",this_season_foa_run,"_",this_season_pt,"_FlameLengths/",
                                     this_season_foa_run, "_", this_season_pt, "_", this_season_scen, "_FlameLengths_FireID_",
                                     this_season_fireIDs, ".tif")
  print(paste0("Reading in Flame Length tifs for Season ", each_season,"..."))
  this_season_FL_stack <- rast(lapply(this_season_FL_filenames, align_and_stack_tifs))
  
  #If there is one or fewer fires in the season, use the process_single_fire_season function 
  if(length(this_season_FL_filenames) <= 1){
    process_single_fire_season(this_season_AD_stack, this_season_FL_stack, each_season, wd, this_season_fireIDs)
  } else { #Otherwise, read in the perimeters sqlite database and fetch this season's fire perimeters
    con <- dbConnect(SQLite(), dbname = paste0(wd,"/",foa_run, "_", this_season_pt[1], "_Perimeters.sqlite"))
    query1 <- paste("SELECT * FROM perimeters WHERE fire_id IN (", toString(this_season_fireIDs),")")
    season_fire_perims <- dbGetQuery(con, query1)
    ref_sys <- dbGetQuery(con, "SELECT * FROM spatial_ref_sys")
    dbDisconnect(con)
    #Convert these perimeters to an sf object
    season_fire_perims$GEOMETRY <- sf::st_as_sfc(structure(as.list(season_fire_perims$GEOMETRY),class="blob"), crs = ref_sys$srtext)
    season_fire_perims <- st_as_sf(season_fire_perims)
    #Plot to confirm they loaded correctly. There are multiple variables, so just plot the first (the object ID).
    plot(season_fire_perims, col = "black", max.plot = 1)
    
    #Check whether any of the fire perimeters for this season overlap. The result is a matrix of logical outcomes.
    overlap_matrix <- st_intersects(season_fire_perims, sparse = FALSE)
    #Use the function find_overlap_indices to determine which fires overlap with which other fires.
    overlap_indices <- find_overlap_indices(overlap_matrix)
    print(paste0("There are ",length(overlap_indices), " cases where fires overlap in season ", each_season, "."))
    
    #If no fires overlap, use the merge_and_write_rasters function to export the 3-band season tif.
    if(length(overlap_indices) == 0){
      merge_and_write_rasters(this_season_AD_stack, this_season_FL_stack, each_season, wd, this_season_fireIDs)
    } else if(length(overlap_indices) >= 1){ #If there's at least one case of overlap, use the function process_overlapping_fires.
      process_overlapping_fires(this_season_AD_stack, this_season_FL_stack, each_season, wd, overlap_indices, season_fire_perims, ref_sys, firelists, foa_run, align_raster, mask, app, extract, global, dbConnect, dbGetQuery, dbDisconnect)
    }
  } 
}


#STEP 5: Process seasons in parallel across multiple cores
#############################################################

# Define the processing function to be run in parallel
process_single_season <- function(each_season, firelists, wd, foa_run, align_and_stack_tifs, dbConnect, dbGetQuery, dbDisconnect, merge_non_na_first, count_non_na, app, global, mask, extract, align_raster, writeRaster) {
  process_fire_season(each_season, firelists, wd, foa_run, align_and_stack_tifs, dbConnect, dbGetQuery, dbDisconnect, merge_non_na_first, count_non_na, app, global, mask, extract, align_raster, writeRaster)
}

# Set up parallel processing plan (adjust workers as needed)
plan(multicore) # For local parallelism; use plan(cluster, workers = num_workers) for a cluster setup

# Get unique seasons
unique_seasons <- unique(firelists$Season)

# Process each season in parallel
results <- future_map(
  unique_seasons,
  process_single_season,
  each_season = each_season,
  firelists = firelists,
  wd = wd,
  foa_run = foa_run
)

# Combine all results
final_results <- do.call(rbind, results)

