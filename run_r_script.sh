#!/bin/bash
# Note: Before running this on Linux, you'll need to run this line to remove stupid Windows characters:
# sed -i -e 's/\r$//' rename_tifs.sh

#SBATCH --job-name=contifs
#SBATCH --partition=math-alderaan
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --output=FOA2c_r10_comb_AD_tifs_slurmlog_%A.out
#SBATCH --error=FOA2c_r10_comb_AD_tifs_error_%A.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=laurel.sindewald@ucdenver.edu

# Change to the directory containing your R script
cd /data001/projects/sindewal/FSim_post_processing/
# Make a directory to hold the merged tifs
mkdir ../okwen_foa2c_r10/SeasonFires_merged_tifs/

# Activate conda environment 
source /data001/projects/sindewal/anaconda3/bin/activate r_env

# Run the R script
/data001/projects/sindewal/anaconda3/envs/r_env/bin/Rscript FSim_post_processing1_merge_fire_tifs_v10_refact_slurm.R --foa_lcp_path ./_inputs/lcp/FOA2c_LCG_LF2022_FBFM40_230_120m.tif --working_directory /data001/projects/sindewal/okwen_foa2c_r10/ --foa_run FOA2c_r10 --scenario "" --run_timepoint baseline_time0 --number_of_seasons 20000 --seasons_in_part 2500 --number_of_parts 8 --first_season 1 --last_season 2500 
