#!/bin/bash

#blast to protein fasta
for filename in ../zipped_fasta/*.fna; do
    tblastn -query sna_protein.fasta -subject "$filename" -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -out $(basename "$filename" .fna)_tblastn
done

for filename in ./*_tblastn; do
    # Sort by the 11th column, numerically
    sort -k 11,11 -g "$filename" > "${filename}_sorted"

    # Extract the top line to represent the best hit
    head -n 2 "${filename}_sorted" > "${filename}_sorted_subset"

    # Move the best hit to a final output file and clean up intermediate files
    mv "${filename}_sorted_subset" "${filename}_best_hit"
    rm "${filename}_sorted"
    rm "$filename"
done

for filename in ./*_best_hit; do
    awk '{
    OFS="\t";
    strand = ($9 <= $10) ? "+" : "-";
    start = ($9 <= $10) ? $9 : $10;
    end = ($9 <= $10) ? $10 : $9;
    print $2, "tblastn", "similarity", start, end, $12, strand, ".", "ID="$1";evalue="$11";bitscore="$12
    }' "$filename" > "$filename".gff
    rm "$filename"
done

#find where the best hit intersect
export PATH=/programs/bedtools2-2.29.2/bin:$PATH
for filename in ./*_best_hit.gff; do
    bedtools intersect -a ../gff_files/$(basename "$filename" _tblastn_best_hit.gff).gtf -b "$filename" -wa > $(basename "$filename" _best_hit.gff)_overlap
    rm "$filename"
done

#isolate all the mRNA strands

export PATH=/programs/bedops-2.4.35/bin:$PATH
for filename in ./*_overlap; do
    grep -w "mRNA" "$filename" > "$filename"_mRNA
    rm "$filename"
done

# find where mRNA strands are and modify annotations
# in actually pipeline, remove all unnecessary files

#!/bin/bash

for filename in ./*_mRNA; do
    # Define the input files
    file_a=../gff_files/$(basename "$filename" _tblastn_overlap_mRNA).gtf # replace with your actual file name
    echo "$file_a"
    file_b="$filename"  # replace with your actual file name
    echo "$file_b"

    # Create a temporary file to store results
    temp_file=$(mktemp)

    # Use awk to process both files, passing temp_file as a variable
    awk -v temp_file="$temp_file" 'BEGIN { OFS="\t" }  # Set output field separator to tab
        # Step 1: Read all lines of file_b and store the 9th column values in an array
        NR == FNR { 
            match_column_b[$9]; 
            next 
        }

        # Step 2: For lines in file_a, check if the 9th column exists in file_b
        {
            # If the 9th column matches, add the suffix ";name=real-sna;"
            if ($9 in match_column_b) {
                $9 = $9";name=real-sna;"
            }

            # Print the modified or original line to the temporary file
            print > temp_file
        }
    ' "$file_b" "$file_a"

    # Replace file_a with the modified content
    mv "$temp_file" "$file_a"
    rm "$filename"

done




for filename in ../gff_files/*.gtf; do
    output_file="bed_genes/$(basename "$filename" .gtf)_genes.gtf"
    > "$output_file"
    echo "$filename"
    
    awk '
    BEGIN { in_gene = 0; strand = ""; cds_line = ""; last_cds_line = "" }
    {
        if ($3 == "mRNA") {
            # Output the last CDS line if previous mRNA was on the - strand
            if (in_gene && strand == "-" && last_cds_line != "") {
                print last_cds_line > output_file
            }

            # Reset variables for a new mRNA block
            in_gene = 0
            cds_line = ""
            last_cds_line = ""
            strand = ""

            # Check if this mRNA is "real-sna"
            if (index($0, "real-sna") != 0) {
                print "sna" > "/dev/stderr"
                in_gene = 1
                strand = $7
            }
        }
        
        # If inside a "real-sna" mRNA annotation and the line is a CDS
        else if (in_gene && $3 == "CDS") {
            if (strand == "+") {
                # For + strand, output the first CDS line found
                if (cds_line == "") {
                    cds_line = $0
                    print cds_line > output_file
                }
            } else if (strand == "-") {
                # For - strand, store the last CDS line encountered
                last_cds_line = $0
            }
        }
    }
    END {
        # After the final mRNA block, append the last CDS line if strand is -
        if (in_gene && strand == "-" && last_cds_line != "") {
            print last_cds_line > output_file
        }
    }
    ' output_file="$output_file" "$filename"
done

export PATH=/programs/bedops-2.4.35/bin:$PATH
for filename in bed_genes/*.gtf; do
    gff2bed < "$filename" > bed_genes/$(basename "$filename" .gtf).bed
    rm "$filename"
done

for filename in bed_genes/*.bed; do
    awk '!seen[$2, $3]++' "$filename" > "$filename"_temp_file && mv "$filename"_temp_file "$filename"
done

for filename in ../zipped_fasta/*.fna; do
    faidx "$filename" -i chromsizes -o fna_sizes/$(basename "$filename" .fna).fna_sizes
done

for filename in bed_genes/*.bed; do
    awk 'BEGIN { OFS = "\t" }
    {
        if ($6 == "-") {
            $2 = $3
            $3 = $2 + 500
        } else if ($6 == "+") {
            $3 = $2
            $2 = $3 - 500
            if ($2 < 0) {
            $2 = 0
            }
        }
        print $0
    }' "$filename" > 500_upstream/minus/$(basename "$filename" .bed)_minus_500_upstream

    grep '+' 500_upstream/minus/$(basename "$filename" .bed)_minus_500_upstream > 500_upstream/plus/$(basename "$filename" .bed)_plus_500_upstream
    sed -i '/+/d' 500_upstream/minus/$(basename "$filename" .bed)_minus_500_upstream
done


for filename in 500_upstream/minus/*_upstream; do
    echo "$filename"
#get fasta for 500 bp windows
    bedtools getfasta -fo "$filename".fasta -fi ../zipped_fasta/$(basename "$filename" _genes_minus_500_upstream).fna -bed "$filename" -nameOnly
    mv "$filename".fasta 500_upstream/minus/fasta
#find motif
    export PATH=/programs/seqkit-0.15.0:$PATH
    seqkit locate --ignore-case -p "CAGGTAG" -p "tAGGTAG" -p "CAGGcAG" 500_upstream/minus/fasta/$(basename "$filename").fasta >"$filename"_canonical_tagteam 
    mv "$filename"_canonical_tagteam 500_upstream/minus/canonical_tagteam
#remove duplicates in file
    awk '!seen[$2, $3]++' 500_upstream/minus/canonical_tagteam/$(basename "$filename")_canonical_tagteam > 500_upstream/minus/canonical_tagteam/$(basename "$filename")_canonical_tagteam_fin
#turn file into .csv format
    sed 's/ \+/,/g' 500_upstream/minus/canonical_tagteam/$(basename "$filename")_canonical_tagteam_fin > 500_upstream/minus/canonical_tagteam/csv/$(basename "$filename").csv
done

for filename in 500_upstream/plus/*_upstream; do
    echo "$filename"
#get fasta for 500 bp windows
    bedtools getfasta -fo "$filename".fasta -fi ../zipped_fasta/$(basename "$filename" _genes_plus_500_upstream).fna -bed "$filename" -nameOnly
    mv "$filename".fasta 500_upstream/plus/fasta
#find motif
    export PATH=/programs/seqkit-0.15.0:$PATH
    seqkit locate --ignore-case -p "CAGGTAG" -p "tAGGTAG" -p "CAGGcAG" 500_upstream/plus/fasta/$(basename "$filename").fasta >"$filename"_canonical_tagteam
    mv "$filename"_canonical_tagteam 500_upstream/plus/canonical_tagteam
#remove duplicates in file
    awk '!seen[$2, $3]++' 500_upstream/plus/canonical_tagteam/$(basename "$filename")_canonical_tagteam > 500_upstream/plus/canonical_tagteam/$(basename "$filename")_canonical_tagteam_fin
#turn file into .csv format
    sed 's/ \+/,/g' 500_upstream/plus/canonical_tagteam/$(basename "$filename")_canonical_tagteam_fin > 500_upstream/plus/canonical_tagteam/csv/$(basename "$filename").csv
done