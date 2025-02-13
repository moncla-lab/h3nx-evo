Downloaded 12/17/2024
Last ran on this download on 12/19/2024

This is a pipeline to go from downloaded ncbi and/or gisaid data to a merged dataset with already existing sequences. This requires you to have a parent folder that has both ncbi and gisaid subfolders with fasta files (and gisaid csv metadata file if using). all other folders are made after running this pipeline

DATA_PULLS FOLDER:
The parent folder where the gisaid and ncbi data is stored. The fasta file in both folders is parsed out into segments (parsed), then each parsed file undergoes QC (QC), which is then deduped (deduped). The last step is to standardize the date across all 8 segments (consistent) in case dates become inconsistent during the deduping step. 

The gisaid and ncbi data is then merged (merged folder), and the same deduping (merged/deduped) and standardizing is performed (merged/consistent).

MASTER_MERGED FOLDER:
This houses the merged fasta files (newly downloaded with already existing)
The master_merged/consistent is what should now be used in new builds