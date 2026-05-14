#!/bin/bash

in_sex_lethal=0
cds_line=""
strand=""

# Loop through all .gtf files
for filename in *.gtf ; do
    # Create the output file for the current GTF file
    output_file="bed_genes/$(basename "$filename" .gtf)_genes.gtf"
    > "$output_file"  # Ensure the output file is empty at the start

    # Remove the temp file if it exists from a previous run
    [ -f temp_cds_last.txt ] && rm temp_cds_last.txt

    echo "Processing file: $filename"  # Debugging message
    
    while IFS= read -r line; do
        # Check if the line contains mRNA information and reset flags
        if [[ "$line" == *"mRNA"* ]]; then
            in_sex_lethal=0  # Reset the flag for each new mRNA
            strand=""  # Reset the strand variable
            
            # Check if the mRNA line corresponds to the sex-lethal product
            if [[ "$line" == *"sex-lethal"* || "$line" == *"Sxl"* ]]; then
                in_sex_lethal=1
                strand=$(echo "$line" | awk '{print $7}')  # Capture strand (+ or -)
                echo "Detected sex-lethal mRNA. Strand: $strand"  # Debugging message
            fi
        fi
    
        # If inside a sex-lethal mRNA annotation, check for CDS lines
        if [[ $in_sex_lethal -eq 1 && "$line" == *"CDS"* ]]; then
            cds_line="$line"
            
            if [[ "$strand" == "+" ]]; then
                # For + strand, write the first CDS line (if no CDS exists yet)
                if [[ -z $(grep "CDS" "$output_file") ]]; then
                    echo "Writing first CDS for + strand: $cds_line"  # Debugging message
                    echo "$cds_line" >> "$output_file"
                fi
            elif [[ "$strand" == "-" ]]; then
                # For - strand, store the last CDS line encountered in temp file
                echo "Storing CDS for - strand: $cds_line"  # Debugging message
                echo "$cds_line" > temp_cds_last.txt
            fi
            # Append the last CDS line for the - strand (if any) to the output file
            if [[ -f temp_cds_last.txt ]]; then
                echo "Appending last CDS line for - strand to $output_file"  # Debugging message
                cat temp_cds_last.txt >> "$output_file"
                rm temp_cds_last.txt
            fi
        fi
    done < "$filename"  
done

#worked
#!/bin/bash

in_sex_lethal=0
cds_line=""
strand=""
last_cds_line=""

# Loop through all .gtf files
for filename in *.gtf ; do
    # Create the output file for the current GTF file
    output_file="bed_genes/$(basename "$filename" .gtf)_genes.gtf"
    > "$output_file"  # Ensure the output file is empty at the start

    while IFS= read -r line; do
        # Check if the line contains mRNA information and reset flags
        if [[ "$line" == *"mRNA"* ]]; then
            # If this is a new mRNA block and the previous block was a - strand, append the last CDS line
            if [[ $in_sex_lethal -eq 1 && "$strand" == "-" && -n "$last_cds_line" ]]; then
                echo "$last_cds_line" >> "$output_file"
            fi

            # Reset for the new mRNA chunk
            in_sex_lethal=0
            last_cds_line=""
            strand=""

            # Check if the mRNA line corresponds to the sex-lethal product
            if [[ "$line" == *"sex-lethal"* || "$line" == *"Sxl"* ]]; then
                in_sex_lethal=1
                strand=$(echo "$line" | awk '{print $7}')  # Capture strand (+ or -)
            fi
        fi
    
        # If inside a sex-lethal mRNA annotation, check for CDS lines
        if [[ $in_sex_lethal -eq 1 && "$line" == *"CDS"* ]]; then
            cds_line="$line"
            
            if [[ "$strand" == "+" ]]; then
                # For + strand, write the first CDS line (if no CDS exists yet)
                if [[ -z $(grep "CDS" "$output_file") ]]; then
                    echo "$cds_line" >> "$output_file"
                fi
            elif [[ "$strand" == "-" ]]; then
                # For - strand, store the last CDS line encountered
                last_cds_line="$cds_line"
            fi
        fi
    done < "$filename"
    
    # After the final mRNA chunk in the file, append the last CDS line for - strand
    if [[ $in_sex_lethal -eq 1 && "$strand" == "-" && -n "$last_cds_line" ]]; then
        echo "$last_cds_line" >> "$output_file"
    fi
done


