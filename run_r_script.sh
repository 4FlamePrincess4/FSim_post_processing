#!/bin/bash
# Note: Before running this on Linux, you'll need to run this line to remove stupid Windows characters:
# sed -i -e 's/\r$//' rename_tifs.sh

#SBATCH --job-name=contifs
#SBATCH --partition=math-alderaan
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=16
#SBATCH --output=FOA2c_r10_tofu_comb_AD_tifs_slurmlog_%A.out
#SBATCH --error=FOA2c_r10_tofu_comb_AD_tifs_error_%A.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=laurel.sindewald@ucdenver.edu

# Load necessary modules 
R

# Define the list of R packages needed
R_PACKAGES=("tidyverse" "sf" "terra" "tidyterra" "RSQLite" "future" "furrr" "log4r")

# Loop through each package and install it if not already installed
for pkg in "${R_PACKAGES[@]}"; do
    library($pkg)
done
echo "Done loading packages."

# Change to the directory containing your R script
cd /data001/projects/sindewal/FSim_post_processing/

# Run the R script
echo "Running R script..."
Rscript FSim_post_processing1_merge_fire_tifs_v10_refactored.R
