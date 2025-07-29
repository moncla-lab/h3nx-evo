#!/bin/bash
#this is a bash script that will clean alignment files that have additional fields separated by a "|"
#will only keep the strain name 

main_folder="path_to_main_folder"


for folder in "$main_folder"/*; do
    if [ -d "$folder" ]; then
        echo "Processing files in folder: $folder"
        
        
        cd "$folder" || exit
        
        
        for fasta_file in *.fa; do
            echo "Processing file: $fasta_file"
            
            python3 /Users/monclalab1/Documents/nonhuman_H3_project/non-human-h3/host_specific/mafft_cleaner.py "$fasta_file" "${fasta_file%.fa}_aligned.fa"
            
            echo "Cleaned alignment created: ${fasta_file%.fa}_aligned.fasta"
        done
        
        
        cd "$main_folder" || exit
    fi
done