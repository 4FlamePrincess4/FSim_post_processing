#!/bin/bash
# Note: Before running this on Linux, you'll need to run this line to remove Windows characters:
# sed -i -e 's/\r$//' rename_tifs_conditional.sh

# Check if the directory argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

# Assign the provided directory argument to a variable
directory="$1"

# Loop through each file in the specified directory
for file in "$directory"/*; do
    # Check if the file is a regular file and exclude directories
    if [ -f "$file" ] && [[ ! "$file" =~ _$ ]]; then
        # Get the base name of the file (without directory path)
        base_file=$(basename "$file")
        # Calculate the length of the filename (including the extension)
        filename_length=${#base_file}
        # Check if the filename (including extension) is less than 26 characters
        if [ "$filename_length" -lt 26 ]; then
            # Get the parent directory name
            current_dir=$(basename "$directory")
            # New filename with parent directory name appended
            new_filename="${current_dir}_$base_file"
            # Rename the file by appending parent directory name and underscore to the beginning
            mv "$file" "$directory/$new_filename"
            # Print the resulting filename
            echo "Renamed '$base_file' to '$new_filename'"
        else
            echo "Skipped '$base_file' (filename length >= 26)"
        fi
    fi
done

echo "Done!"
