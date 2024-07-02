#!/bin/bash
# Note: Before running this on Linux, you'll need to run this line to remove stupid Windows characters:
# sed -i -e 's/\r$//' rename_tifs.sh

#SBATCH --job-name=contifs
#SBATCH --partition=math-alderaan
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --output=FOA2c_r10_tofu_comb_AD_tifs_slurmlog_%A.out
#SBATCH --error=FOA2c_r10_tofu_comb_AD_tifs_error_%A.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=laurel.sindewald@ucdenver.edu

# Change to the directory containing your R script
cd /data001/projects/sindewal/FSim_post_processing/

# Load necessary modules 
R

# Specify the first and last seasons for this job.
first_season <- 1
last_season <- 2500

# Run the R script
echo "Running R script..."
Rscript FSim_post_processing1_merge_fire_tifs_v10_refactored.R
