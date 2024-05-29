library(tidyverse)
library(sf)
library(sp)
library(raster)
library(terra)
library(tidyterra)

#Set the working directory to the specific outputs folder for the run
wd <- setwd("D:/WFSETP/FSim_Outputs/2022_baseline_t0/")

#######################################################################################
# NOTE: To run this code, you need to make sure the following FSim outputs are in the #
#       working directory: the FireSizeList.csv files, the All.tif files, the         #
#       ArrivalDays tifs, the FlameLengths tifs, and the Perimeters.sqlite files.     #
#######################################################################################

run_timepoint <- "baseline_time0"

okwen_foa1c_r7_t0_bp <- raster("./okwen_foa1c_r7/FOA1c_r7_baseline_time0_FullRun_bp.tif")
okwen_foa2c_r10_t0_bp <- raster("./okwen_foa2c_r10/FOA2c_r10_baseline_time0_FullRun_bp.tif")
okwen_foa3c_r3_t0_bp <- raster("./okwen_foa3c_r3/FOA3c_r3_baseline_time0_FullRun_bp.tif")

okwen_extent <- raster("../../Data/Landfire2022_LCGs/Full_OkWen_Landfire2022/Full_OkWen_LCG_LF2022_FBFM40_120m.tif")
okwen_crs <- crs(okwen_extent)
okwen_origin <- origin(okwen_extent)

crs(okwen_foa1c_r7_t0_bp) <- okwen_crs
crs(okwen_foa3c_r3_t0_bp) <- okwen_crs
crs(okwen_foa2c_r10_t0_bp) <- okwen_crs

origin(okwen_foa1c_r7_t0_bp) <- okwen_origin
origin(okwen_foa3c_r3_t0_bp) <- okwen_origin
origin(okwen_foa2c_r10_t0_bp) <- okwen_origin

okwen_foa1c_r7_t0_bp_extended <- raster::extend(okwen_foa1c_r7_t0_bp,
                                                okwen_extent, value=NA)
okwen_foa2c_r10_t0_bp_extended <- raster::extend(okwen_foa2c_r10_t0_bp,
                                                 okwen_extent, value=NA)
okwen_foa3c_r3_t0_bp_extended <- raster::extend(okwen_foa3c_r3_t0_bp,
                                                okwen_extent, value=NA)

okwen_t0_bp_list <- list(okwen_foa1c_r7_t0_bp_extended,okwen_foa2c_r10_t0_bp_extended,okwen_foa3c_r3_t0_bp_extended)
names(okwen_t0_bp_list)[1:2] <- c('x','y')
okwen_t0_bp_list$fun <- sum
okwen_t0_bp_list$na.rm <- TRUE
okwen_t0_bp_list$tolerance <- 4

okwen_t0_bp <- do.call(raster::mosaic, okwen_t0_bp_list)
plot(okwen_t0_bp)
dir.create(paste0("./okwen_", run_timepoint, "_merged_outputs"))
writeRaster(okwen_t0_bp, "./okwen_t0_merged_outputs/okwen_t0_bp.tif", format="GTiff")
