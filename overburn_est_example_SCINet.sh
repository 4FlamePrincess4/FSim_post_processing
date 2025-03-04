#!/bin/bash
# Note: Before running this on Linux, you'll need to run this line to remove stupid Windows characters:
# sed -i -e 's/\r$//' rename_tifs.sh

#SBATCH --job-name=122RO3ovb
#SBATCH --partition=ceres
#SBATCH --time=02-00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --output=/project/wildland_fire_smoke_tradeoff/okwen_foa1c_r16_LF2022_RecOff3_2/FOA1c_r16_2022_RecOff3_2_overburn_summary_slurmlog_%A.out
#SBATCH --error=/project/wildland_fire_smoke_tradeoff/okwen_foa1c_r16_LF2022_RecOff3_2/FOA1c_r16_2022_RecOff3_2_overburn_summary_error_%A.log
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
/home/laurel.sindewald/.conda/envs/r_env2/bin/Rscript FSim_post_processing5_summarize_overburn_removed_Alderaan.R \
--working_directory /project/wildland_fire_smoke_tradeoff/okwen_foa1c_r16_LF2022_RecOff3_2/ \
--season_fires_directory ./SeasonFires_merged_tifs_LF2022_RecOff3_2/ \
--foa_run foa1c_r16 \
--scenario LF2022_RecOff3_2 \
--run_timepoint time2 \
--number_of_seasons 20000 \
--seasons_in_part 5000 \
--number_of_parts 4
