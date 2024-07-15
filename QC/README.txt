MERGED FOLDER:
using the merge() function, the strain and sequences within each segment fasta file from the NCBI and Gisaid data pull is merged with the master fasta file 

DEDUPED FOLDER:
fastaDeDupe() is called here to remove duplicate strains
the longest sequence is kept out of the duplicate pair

CONSISTENT FOLDER:
when deduping, there is the option to also call standardize_dates(). if there is a
duplicate pair, this function will keep the most comprehensive date and apply this date
across each segment. this will ensure the dates are standardized for each strain across 
the genome

Script itself needs more automating/cleanup
