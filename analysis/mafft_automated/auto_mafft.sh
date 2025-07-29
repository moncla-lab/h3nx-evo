#!/bin/bash
#this is a bash script that runs a mafft alignment on a bunch of files at once
#this LACKS the post processing that augur align does (e.g removing insertions relative to a reference)

main_folder="path_to_main_folder"


for folder in "$main_folder"/*; do
    if [ -d "$folder" ]; then
        echo "Processing files in folder: $folder"
        
        #moving into the subfolder
        cd "$folder" || exit
        
        #looping thru each fasta file. make sure the "in *.fa" matches your file extension (.fasta vs .fa)
        for fasta_file in *.fa; do
            echo "Processing file: $fasta_file"
            
            #running mafft on each file; right side is what the output files will be called
            mafft "$fasta_file" > "${fasta_file%.fa}_aligned.fa"
            
            echo "Alignment created: ${fasta_file%.fa}_aligned.fasta"
        done
        
        
        cd "$main_folder" || exit
    fi
done