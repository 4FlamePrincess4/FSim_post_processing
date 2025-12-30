The FSim_post_processing directory contains R and bash scripts to post-process FSim outputs. 
The first number indicates when the script should be launched; some scripts can be launched immediately after FSim finishes, while others are dependent on prior post-processing scripts finishing. 
The second number indexes the scripts in order of importance for each step. 

The R subdirectory contains the following post-processing R scripts.
The bash subdirectory contains scripts to launch the post-processing scripts as slurm jobs in a cluster environment. 
The number of the bash script corresponds to the R script with the same number.

1. FSim_post_processing1_1_merge_run_parts_v2_gdb.R
   This script combines run parts. It 1) merges the FireSizeList.csv files, 2) takes a weighted average of the burn probability and conditional flame length probability rasters (weighted by the seasons in each part), 
   3) reads in the perimeter blobs (binary large objects), merges them into a single vector object, then saves them to an ESRI geodatabase, and 
   4) reads in the ignition points, merges them into a single vector object, then saves them to an ESRI geodatabase.
   --dependency=afterany: FSim runs 

2. FSim_post_processing1_1_2_rename_tifs_conditional.sh
   This script iterates through the ArrivalDays, ArrivalTimes, and FlameLengths fire tif subdirectories and checks the lengths of the file names. If the names are less than or equal to 26 characters in length, 
   the script appends the subdirectory name to the front of each fire tif name. This script is now deprecated, but is currently necessary for scripts 2.2. through 2.7. 
   --dependency=afterany: FSim runs

2.1. FSim_post_processing1_2_merge_fire_tifs_v11_2_Scenario_refact_of_the_Nordgren.R
   This script combines fire tifs (ArrivalDays and FlameLengths) into a 3-band SeasonFire summary raster (Fire IDs, Arrival Days, and Flame Lengths). In the process, the script checks for cases where two fires burn the same area. 
   We assume overburn such as this should not occur, so we use the ArrivalDays raster data to determine which pixels to keep from which fire. In cases where a fire ignites in a previously burned area, we delete the entire fire. 
   This version was made to run in a cluster environment (with options passed via optparse) with a scenario label. The version uses Bryce Nordgren's accumulator method, where we start with blank rasters and add in fire tifs one at a time.
   This versiond does not require fire tif files to be renamed.
   --dependency=afterany: FSim runs 

2.2. FSim_post_processing1_2_merge_fire_tifs_v11_refact_of_the_Nordgren.R
   This version of FSim_post_processing1_2 also runs in a cluster environment and uses Bryce's method, but does not include a scenario label.
   This version currently requires fire tif files to be renamed.
   --depenency=afterok: rename fire tifs

2.3. FSim_post_processing1_2_merge_fire_tifs_v12_Nordgren_verbose_slurm.R
   This version is identical to number 2.2, above, except that it outputs more verbose comments about precisely which fires are kept, corrected, etc. 
   This version currently requires fire tif files to be renamed.
   --depenency=afterok: rename fire tifs

2.4. FSim_post_processing1_2_merge_fire_tifs_v12_Scenario_Nordgren_verbose_slurm.R
   This version is identical to number 2.1, above, except that it outputs more verbose comments about precisely which fires are kept, corrected, etc.
   This version currently requires fire tif files to be renamed.
   --depenency=afterok: rename fire tifs

2.5. FSim_post_processing1_2_merge_fire_tifs_v13_Nordgren_verbose_titan.R
   This version is identical to number 2.3, except that it is designed to run on a stand-alone PC rather than in a cluster environment. 
   This version currently requires fire tif files to be renamed. 
   --depenency=afterok: rename fire tifs

2.6. FSim_post_processing1_2_merge_fire_tifs_v14_Nordgren_fast_titan.R
   This version is identical to number 2.2, except that it is designed to run on a stand-alone PC rather than in a cluster environment.
   This version currently requires fire tif files to be renamed.
   --depenency=afterok: rename fire tifs

2.7. FSim_post_processing1_2_merge_fire_tifs_v15_Nordgren_verbose_titan_onepart.R
   This version assumes that a single FSim run was performed - no parallelization across seasons. This version is otherwise identical to 2.6 except that it is also verbose in comments. 
   This version currently requires fire tif files to be renamed.
   --depenency=afterok: rename fire tifs

3. FSim_post_processing2_1_select_fires_generalized2.R
   This script bootstraps 1000 or 5-year sets, summarizes the distribution of area burned, then selects the 20th, 60th, and 90th decile sets. These sets are later used to simulate fire disturbances in our longitudinal FVS simulations.
   --dependency=afterok: merge run parts

4. FSim_post_processing2_2_summarize_overburn_removed_Alderaan.R
   This script compares acres reported in the merged FireSizeLists to those in the SeasonFire tifs, to estimate overburn removed in each season. The script outputs a text report as well as a csv with summary statistics.
   --dependency=afterok: merge fire tifs

5. FSim_post_processing2_3_recalc_prob_rasters_accumulator_cluster_parallel_5.R
   This script reads in each of the SeasonFire rasters and recalculates the burn probability and conditional flame length probability rasters, accounting for overburn deleted. 
   --dependency=afterok: merge fire tifs

6. FSim_post_processing2_4_mergeFOAs.R
   This script merges the probability rasters across FOAs. The burn probability rasters are summed and the conditional flame length probability rasters are combined with a weighted average (by seasons). 
   The merged conditional flame length probability rasters are multiplied by the merged burn probability raster to output the unconditional flame length probability rasters.
   --dependency=afterok: merge run parts

7. FSim_post_processing2_5_zip_firetif_subdirectories.sh
   This script compresses all directories ending in "ArrivalDays", "ArrivalTimes", "FlameLengths", or "gdb" into zip files.
   --dependency=afterok: merge run parts

7.1. FSim_post_processing2_5_2_auto_compress_tif_folders_conditional.sh
   This script first checks whether all of the fire tif names in subdirectories "ArrivalDays", "ArrivalTimes", and "FlameLengths" are greater than 26 characters, then compresses the subdirectories into zip files.
   This script is now deprecated.

8. FSim_post_processing2_6_delete_firetif_subdirectories.sh
   This script deletes any directories ending in "ArrivalDays", "ArrivalTimes", "FlameLengths", or "gdb" after first checking that there is a corresponding .zip file.
   --dependency=afterok: zip fire tifs
