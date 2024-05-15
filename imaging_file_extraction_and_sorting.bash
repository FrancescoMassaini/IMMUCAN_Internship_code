#!/bin/bash

# Array of cohort prefixes to be processed
cohort_prefixes=("IMMU-BC1" "IMMU-BC2" "IMMU-BC3")  # Add more cohort prefixes as needed
if_types=(IF1 IF2 IF3)  # List of IF types, can be adjusted as needed
search_dir="/home/francesco.massaini/Desktop/IMMUCAN_data/BC2/01_Imaging/IF/cell_properties/"
secondary_pattern="*#_cells_properties2_#*"

# Create a function to handle the movement of files based on the specific patterns
move_files_based_on_pattern() {
    local file=$1
    local cohort_prefix=$2
    local base_file=$(basename "$file" .gz)
    local decompressed_file="${base_file%.gz}"
    local secondary_pattern=$3

    # Extract the file
    echo "Extraction of: $file"
    gzip -dc "$file" > "$search_dir/$decompressed_file"

    # Extract cohort name (e.g., BC1 from IMMU-BC1) for directory structuring
    local cohort_name="${cohort_prefix#IMMU-}"

    # Process each IF type and move files accordingly
    for if_type in "${if_types[@]}"; do
        local target_dir="$search_dir/$cohort_name/$if_type"
        if [[ $decompressed_file == *-${if_type}-* ]]; then
            mkdir -p "$target_dir"
            # Check if the file matches the secondary pattern and determine the correct target sub-directory
            if [[ $decompressed_file == $secondary_pattern ]]; then
                local area_target_dir="$target_dir/Cell_Area"
                mkdir -p "$area_target_dir"
                mv "$search_dir/$decompressed_file" "$area_target_dir"
                echo "Moved $decompressed_file to $area_target_dir"
            else
                mv "$search_dir/$decompressed_file" "$target_dir/"
                echo "Moved $decompressed_file to $target_dir"
            fi
        fi
    done
}

# Process files for each cohort prefix
for cohort_prefix in "${cohort_prefixes[@]}"; do
    echo "Processing files for cohort: $cohort_prefix"
    find "$search_dir" -type f -name "${cohort_prefix}*.tsv.gz" -print0 | while IFS= read -r -d $'\0' file; do
        move_files_based_on_pattern "$file" "$cohort_prefix" "$secondary_pattern"
    done
done
