#!/bin/bash
# Note: Before running this on Linux, you'll need to run this line to remove stupid Windows characters:
# sed -i -e 's/\r$//' rename_tifs.sh

#SBATCH --job-name=120cumbp
#SBATCH --partition=math-alderaan
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --output=FOA1c_r16_2020_cumulative_bp_slurmlog_%A.out
#SBATCH --error=FOA1c_r16_2020_cumulative_bp_error_%A.log
#SBATCH --mail-user=laurel.sindewald@ucdenver.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Change to the directory containing your R script
cd /data001/projects/sindewal/FSim_post_processing/

# Activate conda environment 
source /data001/projects/sindewal/anaconda3/bin/activate r_env2

# Run the R script
/data001/projects/sindewal/anaconda3/envs/r_env2/bin/Rscript cumulative_bp_cluster.R \
--working_directory /data001/projects/sindewal/okwen_foa1c_r16_LF2020/ \
--foa_run foa1c_r16 \
--scenario LF2020 \
--run_timepoint baseline_time0 \
--last_season 20000 \
--first_season 1 \
--num_random_pixels 100000
