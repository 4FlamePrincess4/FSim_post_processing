#!/bin/bash
# Note: Before running this on Linux, you'll need to run this line to remove stupid Windows characters:
# sed -i -e 's/\r$//' rename_tifs.sh

#SBATCH --job-name=F1RO3rec
#SBATCH --account=wildland_fire_smoke_tradeoff
#SBATCH --partition=ceres
#SBATCH --time=04-00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --output=/project/wildland_fire_smoke_tradeoff/okawen_foa1c_r16_LF2020_RO2020_3/_error/foa1c_r16_LF2020_RO2020_3_recalc_prob_rasters_slurmlog_%A.out
#SBATCH --error=/project/wildland_fire_smoke_tradeoff/okawen_foa1c_r16_LF2020_RO2020_3/_error/foa1c_r16_LF2020_RO2020_3_recalc_prob_rasters_error_%A.log
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
/home/laurel.sindewald/.conda/envs/r_env2/bin/Rscript FSim_post_processing6_recalc_prob_rasters_accumulator_cluster_parallel_5.R \
--working_directory /project/wildland_fire_smoke_tradeoff/okawen_foa1c_r16_LF2020_RO2020_3/ \
--foa_lcp_path ./_inputs/lcp/lf2020_220_foa1c_120m.tif \
--foa_run foa1c_r16 \
--scenario LF2020_RO2020_3 \
--run_timepoint time0
