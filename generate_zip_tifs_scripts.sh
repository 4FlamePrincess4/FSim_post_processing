#!/bin/bash

# Usage: ./generate_rename_tifs_scripts.sh <study_area> <foa_run> <scenario> <run_timepoint> <job_name>

if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <study_area> <foa_run> <scenario> <run_timepoint> <job_name>"
    exit 1
fi

study_area=$1
foa_run=$2
scenario=$3
run_timepoint=$4
job_name=$5

output_dir="/project/wildland_fire_smoke_tradeoff/${study_area}_${foa_run}_${scenario}_${run_timepoint}/_error"
script_name="${study_area}_${foa_run}_${scenario}_${run_timepoint}_zip_tifs.sh"

cat > "$script_name" <<EOL
#!/bin/bash

#SBATCH --job-name=${job_name}
#SBATCH --partition=ceres
#SBATCH --time=02-00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --mail-user=laurel.sindewald@usda.gov
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=${output_dir}/${foa_run}_pt${part}_${scenario}_${run_timepoint}_zip_tifs_slurmlog_%A.out
#SBATCH --error=${output_dir}/${foa_run}_pt${part}_${scenario}_${run_timepoint}_zip_tifs_error_%A.log

# Set the base directory to check and minimum allowed filename length
BASE_DIR="/project/wildland_fire_smoke_tradeoff/${study_area}_${foa_run}_${scenario}_${run_timepoint}"   # Replace with the base directory to check
MIN_LENGTH=26                        # Minimum required filename length (to check whether the files were renamed correctly and completely). 
ERROR_FILE="${SLURM_SUBMIT_DIR}/slurm-${SLURM_JOB_ID}.err"  # SLURM error file path

# Flag to determine if compression should occur
compress_needed=true

# Array to store names of directories with filenames below the required length
directories_below_min_length=()

# Find all subdirectories ending in ArrivalDays or FlameLengths
for SUBDIR in "$BASE_DIR"/*{ArrivalDays,FlameLengths}; do
    # Skip if no matching subdirectories are found
    [[ ! -d "$SUBDIR" ]] && continue

    # Check if any filename in the subdirectory is below the minimum length
    SHORT_FILES=$(find "$SUBDIR" -type f -name "*" -printf '%f\n' | awk -v min_length="$MIN_LENGTH" 'length < min_length')

    if [[ -n "$SHORT_FILES" ]]; then
        echo "Not all filenames in $SUBDIR meet the minimum length of $MIN_LENGTH characters."
        directories_below_min_length+=("$SUBDIR")
        compress_needed=false
    else
        echo "All filenames in $SUBDIR meet the minimum length of $MIN_LENGTH characters."
    fi
done

# Log directories with filenames below the required length to the SLURM error file if found
if [[ ${#directories_below_min_length[@]} -gt 0 ]]; then
    echo "Directories with filenames below the required length:" > "$ERROR_FILE"
    for DIR in "${directories_below_min_length[@]}"; do
        echo "$DIR" >> "$ERROR_FILE"
    done
fi

# Only compress directories if all filenames meet the minimum length requirement
if [[ "$compress_needed" == true ]]; then
    echo "All filenames in ArrivalDays and FlameLengths subdirectories meet the minimum length of $MIN_LENGTH characters. Compressing directories..."
    
    # Loop through subdirectories ending in ArrivalDays, FlameLengths, or ArrivalTimes to create individual zip files
    for DIR in "$BASE_DIR"/*{ArrivalDays,FlameLengths,ArrivalTimes}; do
        # Skip if no matching subdirectories are found
        [[ ! -d "$DIR" ]] && continue

        ZIP_NAME="$BASE_DIR/$(basename "$DIR").zip"  # Save the zip file in the BASE_DIR
        zip -r "$ZIP_NAME" "$DIR" 2>/dev/null
        echo "Compressed $DIR into $ZIP_NAME."
    done
else
    echo "Compression skipped due to filenames below the minimum length. Check SLURM error file for details."
fi
EOL

done

echo "Generated SLURM script to compress tifs for ${study_area} ${foa_run}, scenario ${scenario}, run timepoint ${run_timepoint}."
