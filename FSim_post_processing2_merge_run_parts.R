library(tidyverse)
library(sf)
library(sp)
library(raster)
library(terra)
library(tidyterra)
library(RSQLite)
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
#       working directory: the FireSizeList.csv files, the All.tif files, the         #
#       ArrivalDays tifs, the FlameLengths tifs, and the Perimeters.sqlite files.     #
#######################################################################################

#STEP 1: Record run information below 
###############################################
# calculate or parse seasons_per_part
if (is.null(opt$seasons_per_part)) {
  seasons_per_part <- rep(opt$seasons_in_part, opt$number_of_parts)
} else {
  seasons_per_part <- as.integer(unlist(strsplit(opt$seasons_per_part, ",")))
}

#STEP 2: Combine fire lists from all four run parts
###############################################
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
Scenario <- rep(opt$scenario, nrow(firelists))
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

#STEP 3: Merge (average) probability rasters from all parts of each run
###############################################
#Read in the "All" tifs, which are multiband rasters of outputs that include
## 1) burn probability, 2) conditional flame length probability (0-2 ft), 3) conditional flame
## length probability (2-4 ft), 4) CFL prob (4-6 ft), 5) CFL prob (6-8 ft), 6) CFL prob (8-12 ft),
## 7) CFL prob (12+ ft), 8) Conditional Flame Length (ft). 
all_tifs_files <- list.files(path=wd,
                             recursive=F,
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
  plot(mean_p)+title(paste0(opt$foa_run,"_",opt$scenario, "_", opt$run_timepoint,"_FullRun_",layers$layer_name[this_layer]))
  #Write the averaged raster to the working directory
  writeRaster(mean_p,paste0(opt$foa_run,"_",opt$scenario, "_", opt$run_timepoint,"_FullRun_",layers$layer_name[this_layer]),
              format="GTiff", overwrite=TRUE)
  #Re-initialize the p_tifs list
  p_tifs <- list()
}

#STEP 4: Export perimeters for all fires (keep separate by part)
#############################################################
#Read in fire perimeters
perimeter_dbs <- list.files(path=wd,
                            recursive=F,
                            pattern="_Perimeters.sqlite",
                            full.names=T)

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
  st_write(df2,
              paste0("./", "perimeters_", opt$foa_run, "_", opt$scenario, "_", opt$run_timepoint, "_part",i),
              driver= "ESRI Shapefile", append=FALSE)
  #The below options combine the parts into one geojson or shapefile, 
  # but these file types don't support file sizes large enough to store all perimeters.
  # st_write(df2,
  #          paste0("./", "perimeters_", opt$foa_run, "_", opt$scenario, "_", opt$run_timepoint, ".geojson"),
  #          append=TRUE)
  # st_write(df2,
  #          paste0("./", "perimeters_", opt$foa_run, "_", opt$scenario, "_", opt$run_timepoint),
  #          driver= "ESRI Shapefile", append=TRUE)
}
dbDisconnect(con)

#STEP 5: Export ignitions for all fires (combine across parts)
#############################################################
#Read in fire ignitions
point_dbs <- list.files(path=wd,
                        recursive=F,
                        pattern="_Ignitions.sqlite",
                        full.names=T)

i=1

# Initialize an empty vector for accumulating points
pts_vector_all <- NULL

#Loop through each perimeter, manipulate, and then merge all into one file

for (i in seq_along(point_dbs)){
  
  # Extract the run part from the file path (e.g., "pt1" from path)
  run_part <- sub(".*_(pt\\d+)_.*", "\\1", basename(point_dbs[i]))
  
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
  
  # Close the database connection
  dbDisconnect(con)
  
  #create a shapefile from the dataframes ----
  
  #change the x, y geometry fields to lat/long
  colnames(df2)[colnames(df2) == "ignition_x"]<- "lon"
  colnames(df2)[colnames(df2) == "ignition_y"]<- "lat"
  
  # Add run part as a new column in df2
  df2$run_part <- run_part
  
  # Create vector for current points with renamed geometry fields
  pts_vector <- terra::vect(x = df2, geom = c("lon", "lat"), crs = df3$srtext)

  # Combine with accumulated points
  if (is.null(pts_vector_all)) {
    pts_vector_all <- pts_vector
  } else {
    pts_vector_all <- rbind(pts_vector_all, pts_vector)
  }
}
  
#save the file
writeVector(pts_vector_all, paste0("./", "ignitions_", opt$foa_run, "_", opt$scenario, "_", opt$run_timepoint),
              filetype= "ESRI Shapefile")
