#!/bin/bash
#SBATCH --job-name=1ziptifs
#SBATCH --partition=math-alderaan
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --output=FOA1c_r16_2016_ziptifs_slurmlog_%A.out
#SBATCH --error=FOA1c_r16_2016_ziptifs_error_%A.log
#SBATCH --mail-user=laurel.sindewald@ucdenver.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Set the base directory to check and minimum allowed filename length
BASE_DIR="/path/to/base_directory"   # Replace with the base directory to check
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
