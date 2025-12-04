#!/bin/bash
#!/bin/bash
# Batch decompresses NetCDF files with COMPRESSED_SNOW_ prefix to SNOW_ prefix
#
# Usage: bash script.sh /path/to/directory
#
# Check if a path is provided
if [ $# -eq 0 ]; then
    echo "Please provide a directory path"
    echo "Usage: $0 /path/to/directory"
    exit 1
fi

# Get the directory path from the first argument
DIR="$1"

# Ensure the directory exists
if [ ! -d "$DIR" ]; then
    echo "Directory does not exist: $DIR"
    exit 1
fi

# Change to the specified directory
cd "$DIR" || exit

# Loop through all files matching the pattern
for file in COMPRESSED_SNOW_*.nc; do
    # Check if files exist
    if [ -e "$file" ]; then
        # Create the new filename by replacing COMPRESSED_SNOW with SNOW
        new_file="${file/COMPRESSED_SNOW/SNOW}"

        # Decompress the file
        nccopy -d0 -s "$file" "$new_file"

        echo "Decompressed: $file -> $new_file"
    fi
done
