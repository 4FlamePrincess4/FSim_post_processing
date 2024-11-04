#!/bin/bash

#SBATCH --job-name=1rentifs
#SBATCH --partition=math-alderaan
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --output=FOA1c_r16_2016_rename_tifs_slurmlog_%A.out
#SBATCH --error=FOA1c_r16_2016_rename_tifs_error_%A.log
#SBATCH --mail-user=laurel.sindewald@usda.gov
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# This script checks each file name in subdirectories ending with ArrivalDays or FlameLengths,
# and if a file name is shorter than 26 characters, it prepends the subdirectory's base name.

# FSim run directory to process
root_directory="/path/to/your/directory"  

# Redirect all output and errors to the SLURM error log file
exec 1>&2

# Check if the root directory exists
if [ ! -d "$root_directory" ]; then
    echo "Error: Directory '$root_directory' does not exist."
    exit 1
fi

# Find all subdirectories ending in ArrivalDays or FlameLengths
find "$root_directory" -type d \( -name "*ArrivalDays" -o -name "*FlameLengths" \) | while read -r subdir; do
    # Get the base name of the subdirectory
    current_dir=$(basename "$subdir")
    
    # Loop through each file in the current subdirectory
    for file in "$subdir"/*; do
        # Check if the file is a regular file and exclude directories
        if [ -f "$file" ] && [[ ! "$file" =~ _$ ]]; then
            # Get the base name of the file (without directory path)
            base_file=$(basename "$file")
            # Calculate the length of the filename (including the extension)
            filename_length=${#base_file}
            # Check if the filename (including extension) is less than 26 characters
            if [ "$filename_length" -lt 26 ]; then
                # New filename with parent directory name appended
                new_filename="${current_dir}_$base_file"
                # Rename the file by appending parent directory name and underscore to the beginning
                mv "$file" "$subdir/$new_filename"
                # Print the resulting filename
                echo "Renamed '$base_file' to '$new_filename' in '$subdir'"
            else
                echo "Skipped '$base_file' in '$subdir' (filename length >= 26)"
            fi
        fi
    done
done

echo "Done renaming files in ArrivalDays and FlameLengths subdirectories!"
