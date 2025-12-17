#!/bin/bash
# Note: Before running this on Linux, you'll need to run this line to remove stupid Windows characters:
# sed -i -e 's/\r$//' rename_tifs.sh

#SBATCH --job-name=120BOzip
#SBATCH --account=wildland_fire_smoke_tradeoff
#SBATCH --partition=ceres
#SBATCH --time=02-00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --output=/project/wildland_fire_smoke_tradeoff/okawen_foa1c_r17_LF2020_TM_baseline_time0/_error/foa1c_r17_LF2020_TM_baseline_time0_ziptifs_slurmlog_%A.out
#SBATCH --error=/project/wildland_fire_smoke_tradeoff/okawen_foa1c_r17_LF2020_TM_baseline_time0/_error/foa1c_r17_LF2020_TM_baseline_time0_ziptifs_error_%A.log
#SBATCH --mail-user=laurel.sindewald@usda.gov
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

BASE_DIR="/project/wildland_fire_smoke_tradeoff/okawen_foa1c_r17_LF2020_TM_baseline_time0"   # Replace with the base directory 
ERROR_FILE="${SLURM_SUBMIT_DIR}/slurm-${SLURM_JOB_ID}.err"  # SLURM error file path

# Loop through subdirectories ending in ArrivalDays, FlameLengths, or ArrivalTimes to create individual zip files
 for DIR in "$BASE_DIR"/*{ArrivalDays,FlameLengths,ArrivalTimes}; do
    # Skip if no matching subdirectories are found
    [[ ! -d "$DIR" ]] && continue

    ZIP_NAME="$BASE_DIR/$(basename "$DIR").zip"  # Save the zip file in the BASE_DIR
    zip -r "$ZIP_NAME" "$DIR" 2>/dev/null
    echo "Compressed $DIR into $ZIP_NAME."
 done
fi
