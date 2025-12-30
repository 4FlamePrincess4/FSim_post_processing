#!/bin/bash
# Note: Before running this on Linux, you'll need to run this line to remove stupid Windows characters:
# sed -i -e 's/\r$//' rename_tifs.sh

#SBATCH --job-name=116overb
#SBATCH --partition=math-alderaan
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --output=FOA1c_r16_2016_overburn_summary_slurmlog_%A.out
#SBATCH --error=FOA1c_r16_2016_overburn_summary_error_%A.log
#SBATCH --mail-user=laurel.sindewald@ucdenver.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Change to the directory containing your R script
cd /data001/projects/sindewal/FSim_post_processing/

# Activate conda environment 
source /data001/projects/sindewal/anaconda3/bin/activate r_env2

# Run the R script
/data001/projects/sindewal/anaconda3/envs/r_env2/bin/Rscript FSim_post_processing5_summarize_overburn_removed_Alderaan.R \
--working_directory /data001/projects/sindewal/okwen_foa1c_r16_LF2016/ \
--season_fires_directory ./SeasonFires_merged_tifs_LF2016/ \
--foa_run foa1c_r16 \
--scenario LF2016 \
--run_timepoint baseline_time0 \
--number_of_seasons 20000 \
--seasons_in_part 5000 \
--number_of_parts 4 
