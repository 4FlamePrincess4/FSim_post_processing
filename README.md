The FSim post-processing R scripts were created to efficiently 1) combine individual fire tifs from FSim runs into a single tif for each fire season while deleting overburn, 2) merge outputs from FSim runs split into multiple parts and output fire perimeter and ignition point shapefiles from the sqlite database outputs, 3) sub-sample seasons to update a landscape with disturbance in FVS or estimate emissions across a range of fire season severities, 4) merge probability raster outputs from FOAs across a full study area, 5) estimate overburn deleted by step 1, and 6) recalculate the burn probability, conditional flame length probability, and unconditional flame length probability rasters from the SeasonFire rasters (accounting for the overburn area deleted). We may not need the select fires step, so that is ignored in the pipeline below. The code has been retained in case it ends up being useful. 

Additional scripts were created to automate processes on a linux cluster or Titan machine. For example, these early iteration scripts were created to rename fire tifs: rename_tifs.ps1, rename_tifs.sh, rename_tifs2.ps1. The linux version of this script was updated to remove the time-consuming process of copying the script to directories and then running them. It was then updated the Linux version to first check whether files had been renamed, and then either rename the files or not and print a statement to an output file: rename_tifs_conditional.sh.

Finally, two scripts were created that automate two time-consuming processes on a linux cluster: 1) auto_rename_tifs_conditional.sh = renaming all the Arrival Day and Flame Length tifs and 2) auto_compress_tif_folders.sh = once those are renamed, compress all of the ArrivalDays, FlameLengths, and ArrivalTimes subdirectories. These scripts are designed to be submitted as slurm jobs dependent on the previous scripts in the pipeline finishing. 

The full pipeline now is: 
1. Run FSim
2. FSim_post_processing2_merge_run_parts.R -> FSim_post_processing4_merge_FOAs.R
2. auto_rename_tifs_conditional.sh -> auto_compress_tif_folders.sh
3. FSim_post_processing1_merge_fire_tifs_v11_Scenario_refact_of_the_Nordgren.R 
4. FSim_post_processing5_summarize_overburn_removed_Alderaan.R
5. FSim_post_processing6_recalc_prob_rasters_accumulator_cluster_parallel_5.R

The R scripts each require a separate bash script to run the R script as a slurm job. Currently there are the run_r_script.sh and run_r_script2.txt template files. 

There are seven versions of FSim_post_processing1 (the script to merge the Arrival Day and Flame Length tifs by season while correcting overburn): 
1. FSim_post_processing1_merge_fire_tifs_v11_refact_of_the_Nordgren.R is a version for Linux clusters that merges tifs with an accumulator method suggested by Bryce Nordgren. For each pixel, the code determines whether the new arrival day values are less than the existing values and, if they are, they replace the value in an accumulator raster. The raster therefore accumulates the lowest arrival day values for fires across the landscape, and the arrival day raster is used as a mask for the flame length and fire ID rasters, which are then added to flame length and fire ID accumulator rasters. This method has lower memory demands than other approaches tried. This version is faster because it doesn't include a check for how many fires overlap at a given location (if overlaps occur). 
2. FSim_post_processing1_merge_fire_tifs_v11_Scenario_refact_of_the_Nordgren.R is identical to the above, except that the code has been adjusted to accommodate a scenario label in the file names. Previously, the naming convention was foa_run_part_file.out. Now, the naming convention is foa_run_part_scenario_file.out.
3. FSim_post_processing1_merge_fire_tifs_v12_Nordgren_verbose_slurm.R is another version for Linux clusters that matches script #1, above, except that it includes a check for how many fires overlap at a given location. It includes more print statements to describe each overlap case and how the overburn was removed. This version is much more memory intensive than scripts 1 or 2, above.
4. FSim_post_processing1_merge_fire_tifs_v12_Scenario_Nordgren_verbose_slurm.R is identical to #3, above, except that the code has been adjusted to accommodate a scenario label in the file names.
5. FSim_post_processing1_merge_fire_tifs_v13_Nordgren_verbose_titan.R is the same as script #3, except that the parallelization code has been set for a single computer instead of a Linux cluster, and runs like a typical R script.
6. FSim_post_processing1_merge_fire_tifs_v14_Nordgren_fast_titan.R is the same as script #1, except that the parallelization code has been set for a single computer instead of a Linux cluster, and runs like a typical R script. 
7. FSim_post_processing1_merge_fire_tifs_v15_Nordgren_verbose_titan_onepart.R is a simplified version of #5 created for users who have not broken their FSim run into parts (chunks of seasons), but ran the seasons together, sequentially, in typical fashion.

These bash scripts are included as templates/examples: 
1. okwen_foa1c_r16_2016_contifs.sh is a script to launch one of the FSim_post_processing1 scripts on a Linux cluster.
2. okwen_foa1c_r16_2016_overburn.sh is a script to launch the FSim_post_processing5 script for estimating overburn on a Linux cluster.
3. RunFSim_Sing_foa1c_r16_pt1_LF2020.sh is a script to submit running FSim as a slurm job. It launches the the singularity container and runs the equivalent of the usual .bat file.
4. RunFSim_foa1c_r16_pt1_LF2020.sh is the Linux equivalent of the usual FSim bat file, which runs FSim referencing a cmdx file, outputs a version file, and logs progress to a _captainslog.txt file.
