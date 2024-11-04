The FSim post-processing R scripts were created to efficiently 1) combine individual fire tifs from FSim runs into a single tif for each fire season while deleting overburn, 2) merge outputs from FSim runs split into multiple parts and output fire perimeter and ignition point shapefiles from the sqlite database outputs, 3) sub-sample seasons to update a landscape with disturbance in FVS or estimate emissions across a range of fire season severities, 4) merge probability raster outputs from FOAs across a full study area, and 5) estimate overburn deleted by step 1. We may not need the select fires step, so that is ignored in the pipeline below. The code has been retained in case it ends up being useful. 

Additional scripts were created to automate processes on a linux cluster or Titan machine. For example, these early iteration scripts were created to rename fire tifs: rename_tifs.ps1, rename_tifs.sh, rename_tifs2.ps1. The linux version of this script was updated to remove the time-consuming process of copying the script to directories and then running them. It was then updated the Linux version to first check whether files had been renamed, and then either rename the files or not and print a statement to an output file: rename_tifs_conditional.sh.

Finally, two scripts were created that automate two time-consuming processes on a linux cluster: 1) renaming all the Arrival Day and Flame Length tifs and 2) once those are renamed, compress all of the ArrivalDays, FlameLengths, and ArrivalTimes subdirectories. These scripts are designed to be submitted as slurm jobs dependent on the previous scripts in the pipeline finishing. 

The full pipeline now is: 
1. Run FSim
2. FSim_post_processing2_merge_run_parts.R -> FSim_post_processing4_merge_FOAs.R
2. auto_rename_tifs_conditional.sh -> auto_compress_tif_folders.sh
3. FSim_post_processing1_merge_fire_tifs_v11_Scenario_refact_of_the_Nordgren.R 
4. FSim_post_processing5_summarize_overburn_removed_Alderaan.R

The R scripts each require a separate bash script to run the R script as a slurm job. Currently there are the run_r_script.sh and run_r_script2.txt template files. Soon, I will add template slurm job scripts for each type of R job since there are fine variations in the slurm scripts depending on whether there are options, etc. 
