#!/bin/bash

set -euo pipefail

# --- Configuration ---
gene_name="sxl" # replace with your gene of interest
protein_fasta="${gene_name}_protein.fasta" # protein sequence fasta file for gene of interest
fasta_dir="../zipped_fasta" #location of fasta files for all species
gtf_dir="../gtf_files" #location of gtf files for all species


mkdir -p ./{bed_genes,500_upstream/{plus/{fasta,canonical_tagteam},minus/{fasta,canonical_tagteam}},results/{plus,minus}}
# The directory is as follows in which this is the file main.sh in the workding directory gene-name
#   gene_name/
#      500_upstream/
#          plus/
#              fasta/
#              canonical_tagteam/
#                csv/ #where final coordinates are located
#           minus/
#              fasta/
#              canonical_tagteam/
#                csv/ #where final coordinates are located
#      bed_genes/
#      results/
#           plus/
#           minus/
#

export PATH=/programs/bedtools2-2.29.2/bin:$PATH
export PATH=/programs/bedops-2.4.35/bin:$PATH
export PATH=/programs/seqkit-0.15.0:$PATH


# --- Step 1: Run tblastn against all fna files ---
# Description: Blasts each species genome to the gene of interest and finds the best hit. Determines where in the genome the best hit corresponds to.
export protein_fasta fasta_dir gtf_dir N

process_file() {
    fasta_file="$1"
    base=$(basename "$fasta_file" .fna)

    tblastn -query "$protein_fasta" -subject "$fasta_file" -outfmt 6 -evalue 1e-5 -max_target_seqs 10 | # blast to the protein sequence
        head -n 1 | # get the best hit located in the first row
        awk '{
            OFS="\t";
            strand = ($9 <= $10) ? "+" : "-";
            start = ($9 <= $10) ? $9 : $10;
            end = ($9 <= $10) ? $10 : $9;
            print $2, "tblastn", "similarity", start, end, $12, strand, ".", "ID="$1";evalue="$11";bitscore="$12
        }' > "${base}_bh.gff"

    bedtools intersect -a "$gtf_dir/${base}.gtf" -b "${base}_bh.gff" -wa > "${base}_bh_overlap" # find where the best hit overlaps in the genome

    rm "${base}_bh.gff"
}

export -f process_file

echo "Running tblastn searches..."
parallel -j "$N" process_file ::: "$fasta_dir"/*.fna

# --- Step 2: Isolate only the mRNA strands in which the best hit overlaps with  ---
echo "Isolating mRNA strands..."
for overlap_file in ./*_overlap; do
    grep -w "mRNA" "$overlap_file" > "${overlap_file}_mRNA"
    rm "$overlap_file"  # Remove the original overlap file after processing
done

# --- Step 3: Modify gtf files to label the best hit gene ---
echo "Modifying annotations..."
# Modify annotations in the mRNA files
for mRNA_file in ./*_mRNA; do
    file_a=""$gtf_dir"/$(basename "$mRNA_file" _bh_overlap_mRNA).gtf"
    file_b="$mRNA_file"
    
# Create a temporary file for modified output
    temp_file=$(mktemp)

    # Use awk to process both files, passing temp_file as a variable
    awk -v gene_name="$gene_name" -v temp_file="$temp_file" 'BEGIN { OFS="\t" }  # Set output field separator to tab
        # Step 1: Read all lines of file_b and store the 9th column values in an array
        NR == FNR { 
            match_column_b[$9]; 
            next 
        }

        # Step 2: For lines in file_a, check if the 9th column exists in file_b
        {
#	If	the	9th	column	matches,	add	the	suffix suffix ";name=real-gene-name;"
            if ($9 in match_column_b) {
            $9	=	$9";name=" gene_name ";"						
            }

            # Print the modified or original line to the temporary file
            print > temp_file
        }
    ' "$file_b" "$file_a"

    # Replace file_a with the modified content
    mv "$temp_file" "$file_a"
    rm "$mRNA_file"

done

# --- Step 4: Find the start codon for the gene  ---
for filename in "$gtf_dir"/*.gtf; do
    output_file="bed_genes/$(basename "$filename" .gtf)_genes.gtf"
    > "$output_file"
    echo "$filename"
    
    awk -v gene_name="$gene_name" '
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

#	Check	if	this	mRNA	is	"real-gene-name"		
if	(index($0,	"gene_name")	!=	0)	{			
                print "gene_name" > "/dev/stderr"
                in_gene = 1
                strand = $7
            }
        }
        
#	If	inside	a	"real-gene-name"	mRNA	annotation	and	the the line is a CDS
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

# --- Step 5: Isolate 500 bp upstream of start codon  ---
for gtf_file in bed_genes/*.gtf; do
    base=$(basename "$gtf_file" .gtf)
    
    # Convert GTF to BED, remove duplicates, and extract 500bp upstream regions
    gff2bed < "$gtf_file" | awk -v base="$base" '
        BEGIN { OFS = "\t" }
        !seen[$2, $3]++ {
            if ($6 == "-") {
                start = $3
                end = start + 500
                print $1, start, end, $4, $5, $6 > "500_upstream/minus/" base "_minus_500_upstream"
            } else if ($6 == "+") {
                end = $2
                start = end - 500
                if (start < 0) start = 0
                print $1, start, end, $4, $5, $6 > "500_upstream/plus/" base "_plus_500_upstream"
            }
        }
    '

    # Remove the original GTF file
    rm "$gtf_file"
done

# --- Step 6: Get the coordinates for the motif locations in base pair region ---
# Define function to process strand
process_strand() {
    strand=$1  # "plus" or "minus"
    input_dir="500_upstream/$strand"

    for filepath in "$input_dir"/*_upstream; do
        echo "Processing $filepath"

        base=$(basename "$filepath")
        fasta_file="$input_dir/fasta/${base}.fasta"
        genome_fasta="$fasta_dir/$(basename "$filepath" _genes_${strand}_500_upstream).fna"
        tagteam_output="$input_dir/canonical_tagteam/${base}_canonical_tagteam"
        final_output="${tagteam_output}_fin"
        csv_output="./results/$strand/${base}.csv"

        # Get FASTA
        bedtools getfasta -fo "$fasta_file" -fi "$genome_fasta" -bed "$filepath" -nameOnly

        # Locate motifs
        seqkit locate --ignore-case -p "CAGGTAG" -p "tAGGTAG" -p "CAGGcAG" "$fasta_file" > "$tagteam_output"

        # Remove duplicate coordinates
        awk '!seen[$2, $3]++' "$tagteam_output" > "$final_output"

        # Convert to CSV format
        sed 's/ \+/,/g' "$final_output" > "$csv_output"
    done
}


export -f process_strand
# Run both strands
parallel process_strand ::: "minus" "plus"

for filename in "$gtf_dir"/*; do
    input_file="$filename"
    output_file="$filename".2

    awk -v gene_name="$gene_name" 'BEGIN { OFS="\t" }
    {
        # Check if "gene_name" is in any column
        for (i=1; i<=NF; i++) {
            if ($i ~ /gene_name/) {
                # Combine all columns from 9 onwards into column 9
                for (j=9; j<=NF; j++) {
                    $9 = $9 " " $j;
                }
                # Set NF to 9, effectively removing all columns after the 9th
                NF = 9;
                break;  # Exit the loop after modifying the columns
            }
        }
        print $0
    }' "$input_file" > "$output_file"
    mv "$output_file" "$input_file"
done

