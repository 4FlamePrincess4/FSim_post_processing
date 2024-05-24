#!/bin/bash
# Note: Before running this on Linux, you'll need to run this line to remove stupid Windows characters:
# sed -i -e 's/\r$//' rename_tifs.sh

#SBATCH --job-name=contifs
#SBATCH --partition=math-alderaan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --output=FOA2c_r10_tofu_comb_AD_tifs_slurmlog_%A.out
#SBATCH --error=FOA2c_r10_tofu_comb_AD_tifs_error_%A.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=laurel.sindewald@ucdenver.edu

# Load necessary modules (if needed)
module load R

# Define the list of R packages needed
R_PACKAGES=("tidyverse" "sf" "sp" "raster" "terra" "tidyterra" "RSQLite")

# Loop through each package and install it if not already installed
for pkg in "${R_PACKAGES[@]}"; do
    if ! R -e "library($pkg)" &>/dev/null; then
        echo "Installing $pkg..."
        R -e "install.packages('$pkg', repos='https://cran.rstudio.com/')"
    else
        echo "$pkg is already installed."
    fi
done
echo "Done loading and installing packages."

# Change to the directory containing your R script
cd /data001/projects/sindewal/FSim_post_processing/

# Run the R script
echo "Running R script..."
Rscript FSim_post-processing1_merge_fire_tifs.R
