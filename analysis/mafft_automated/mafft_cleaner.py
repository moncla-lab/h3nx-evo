
#this is a script that will clean alignment files that have additional fields separated by a "|"
#will only keep the strain name 
#this gets called in the bash script "clean_align.sh"

import sys
import pandas as pd

def mafft_clean(input_file, output_file):

    fasta_data = []

    with open(input_file) as f:
        header = ""
        sequence = ""
        for line in f:
            if line.startswith(">"):
                if header != "":
                    fasta_data.append({"header": header, "sequence": sequence})
                header = line.strip() 
                sequence = ""
            else:
                sequence += line.strip()
        fasta_data.append({"header": header, "sequence": sequence}) #last line

    df = pd.DataFrame(fasta_data)
    df['header'] = df['header'].str.split("|").str[0]

    with open(output_file, "w") as f:
        for index, row in df.iterrows():
            f.write(f"{row['header']}\n")
            f.write(f"{row['sequence']}\n")

if __name__ == "__main__":
	input_file = sys.argv[1]
	output_file = sys.argv[2]
	mafft_clean(input_file, output_file)