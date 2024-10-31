library(tidyverse)
library(sf)
library(sp)
library(raster)
library(terra)
library(tidyterra)

#Set the working directory to the specific outputs folder for the run
wd <- setwd("D:/WFSETP/FSim_Outputs/TTN_LF2020/")

#######################################################################################
# NOTE: To run this code, you need to make sure the following FSim outputs are in the #
#       working directory: the FireSizeList.csv files, the All.tif files, the         #
#       ArrivalDays tifs, the FlameLengths tifs, and the Perimeters.sqlite files.     #
#######################################################################################

run_timepoint <- "baseline_time0"
scenario <- "LF2020"

okwen_foa1c_r16_t0_bp <- raster("./okwen_foa1c_r16/foa1c_r16_LF2020_baseline_time0_FullRun_bp.tif")
okwen_foa2d_r5_t0_bp <- raster("./okwen_foa2d_r5/foa2d_r5_LF2020_baseline_time0_FullRun_bp.tif")
okwen_foa3d_r8_t0_bp <- raster("./okwen_foa3d_r8/foa3d_r8_LF2020_baseline_time0_FullRun_bp.tif")
colville_foa4d_r7_t0_bp <- raster("./colville_foa4d_r7/foa4d_r7_LF2020_baseline_time0_FullRun_bp.tif")
colville_foa5d_r10_t0_bp <- raster("./colville_foa5d_r10/foa5d_r10_LF2020_baseline_time0_FullRun_bp.tif")

okwen_extent <- raster("../../Data/LF2020_LCP_2023_OKAWEN_Colville/LF2022_230_OKAWEN_Colville_LCP-001.tif")
okwen_crs <- crs(okwen_extent)
okwen_origin <- origin(okwen_extent)

crs(okwen_foa1c_r16_t0_bp) <- okwen_crs
crs(okwen_foa2d_r5_t0_bp) <- okwen_crs
crs(okwen_foa3d_r8_t0_bp) <- okwen_crs
crs(colville_foa4d_r7_t0_bp) <- okwen_crs
crs(colville_foa5d_r10_t0_bp) <- okwen_crs

origin(okwen_foa1c_r16_t0_bp) <- okwen_origin
origin(okwen_foa2d_r5_t0_bp) <- okwen_origin
origin(okwen_foa3d_r8_t0_bp) <- okwen_origin
origin(colville_foa4d_r7_t0_bp) <- okwen_origin
origin(colville_foa5d_r10_t0_bp) <- okwen_origin

okwen_foa1c_r16_t0_bp_extended <- raster::extend(okwen_foa1c_r16_t0_bp,
                                                okwen_extent, value=NA)
okwen_foa2d_r5_t0_bp_extended <- raster::extend(okwen_foa2d_r5_t0_bp,
                                                 okwen_extent, value=NA)
okwen_foa3d_r8_t0_bp_extended <- raster::extend(okwen_foa3d_r8_t0_bp,
                                                okwen_extent, value=NA)
colville_foa4d_r7_t0_bp_extended <- raster::extend(colville_foa4d_r7_t0_bp,
                                                okwen_extent, value=NA)
colville_foa5d_r10_t0_bp_extended <- raster::extend(colville_foa5d_r10_t0_bp,
                                                okwen_extent, value=NA)

okwen_t0_bp_list <- list(okwen_foa1c_r16_t0_bp_extended,okwen_foa2d_r5_t0_bp_extended,okwen_foa3d_r8_t0_bp_extended,
                         colville_foa4d_r7_t0_bp_extended, colville_foa5d_r10_t0_bp_extended)
names(okwen_t0_bp_list)[1:2] <- c('x','y')
okwen_t0_bp_list$fun <- sum
okwen_t0_bp_list$na.rm <- TRUE
okwen_t0_bp_list$tolerance <- 4

okwen_t0_bp <- do.call(raster::mosaic, okwen_t0_bp_list)
plot(okwen_t0_bp)
merged_dir <- paste0("./okwen_", run_timepoint, "_", scenario, "_merged_outputs")
dir.create(merged_dir)
writeRaster(okwen_t0_bp, 
            paste0(merged_dir, "/okwen_", run_timepoint, "_", scenario,"merged_bp.tif"), format="GTiff")
