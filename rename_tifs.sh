#!/bin/bash
# Note: Before running this on Linux, you'll need to run this line to remove stupid Windows characters:
# sed -i -e 's/\r$//' rename_tifs.sh

# Loop through each file in the directory
for file in *; do
    # Check if the file is a regular file and exclude directories ending with "_"
    if [ -f "$file" ] && [[ ! "$file" =~ _$ ]]; then
        # Get the parent directory name
        current_dir=$(basename "$(pwd)")
        # New filename with parent directory name appended
        new_filename="${current_dir}_$file"
        # Rename the file by appending parent directory name and underscore to the beginning
        mv "$file" "$new_filename"
        # Print the resulting filename
        echo "Renamed '$file' to '$new_filename'"
    fi
done

echo "Done!"
