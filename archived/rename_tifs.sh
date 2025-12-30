#!/bin/bash
# Note: Before running this on Linux, you'll need to run this line to remove Windows characters:
# sed -i -e 's/\r$//' rename_tifs.sh

# Check if the directory argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

# Assign the provided directory argument to a variable
directory="$1"

# Loop through each file in the specified directory
for file in "$directory"/*; do
    # Check if the file is a regular file and exclude directories ending with "_"
    if [ -f "$file" ] && [[ ! "$file" =~ _$ ]]; then
        # Get the parent directory name
        current_dir=$(basename "$directory")
        # Get the base name of the file (without directory path)
        base_file=$(basename "$file")
        # New filename with parent directory name appended
        new_filename="${current_dir}_$base_file"
        # Rename the file by appending parent directory name and underscore to the beginning
        mv "$file" "$directory/$new_filename"
        # Print the resulting filename
        echo "Renamed '$file' to '$new_filename'"
    fi
done

echo "Done!"
