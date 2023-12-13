#!/bin/bash

# Input GFF file path
input_file="sRNAs_combined.gff"

# Sort the GFF file by chromosome and start position
sorted_file=$(mktemp)
awk '{print $0, $5 - $4}' "$input_file" | sort -k1,1 -k4,4n -k10,10n > "$sorted_file"

# Process the sorted GFF file
output_file="filtered_sRNAs.gff"
prev_line=""
prev_length=0

while IFS= read -r line; do
    curr_chrom=$(echo "$line" | awk '{print $1}')
    curr_start=$(echo "$line" | awk '{print $4}')
    curr_end=$(echo "$line" | awk '{print $5}')
    curr_length=$(echo "$line" | awk '{print $NF}')

    # Compare with the previous line
    if [[ "$prev_line" != "" && "$curr_chrom" == "$(echo "$prev_line" | awk '{print $1}')" ]]; then
        prev_end=$(echo "$prev_line" | awk '{print $5}')
        
        # Check for overlap
        if [[ "$curr_start" -le "$prev_end" ]]; then
            # Overlapping sRNAs, keep the longer one
            if [[ "$curr_length" -gt "$prev_length" ]]; then
                # Replace the previous line with the current line
                prev_line="$line"
                prev_length="$curr_length"
            fi
        else
            # No overlap, write the previous line to the output
            echo "$prev_line" >> "$output_file"
            prev_line="$line"
            prev_length="$curr_length"
        fi
    else
        # First line or different chromosome, update previous line
        prev_line="$line"
        prev_length="$curr_length"
    fi
done < "$sorted_file"

# Write the last line to the output
echo "$prev_line" >> "$output_file"

# Clean up temporary files
rm "$sorted_file"

echo "Filtered sRNAs written to $output_file"
