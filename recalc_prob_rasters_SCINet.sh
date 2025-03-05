#!/bin/bash
# Note: Before running this on Linux, you'll need to run this line to remove stupid Windows characters:
# sed -i -e 's/\r$//' rename_tifs.sh

#SBATCH --job-name=122RO3recalc
#SBATCH --partition=ceres
#SBATCH --time=02-00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --output=/project/wildland_fire_smoke_tradeoff/okwen_foa1c_r16_LF2022_RecOff3_2/FOA1c_r16_2022_RecOff3_2_recalc_prob_rasters_slurmlog_%A.out
#SBATCH --error=/project/wildland_fire_smoke_tradeoff/okwen_foa1c_r16_LF2022_RecOff3_2/FOA1c_r16_2022_RecOff3_2_recalc_prob_rasters_error_%A.log
#SBATCH --mail-user=laurel.sindewald@usda.gov
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Change to the directory containing your R script
cd /project/wildland_fire_smoke_tradeoff/FSim_post_processing/

# Activate conda environment
module load miniconda
source activate r_env2

# Run the R script
/home/laurel.sindewald/.conda/envs/r_env2/bin/Rscript recalc_prob_rasters_cluster.R \
--working_directory /data001/projects/sindewal/okwen_foa1c_r16_LF2020/ \
--foa_run foa1c_r16 \
--scenario LF2020 \
--run_timepoint baseline_time0
