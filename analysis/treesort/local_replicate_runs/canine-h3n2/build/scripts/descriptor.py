import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--descriptor', type=str, required=True, help='path to descriptor.csv')

args = parser.parse_args()

list_of_genes = ["ha", "pb2","pb1","na","np","pa","ns","mp"]

descriptor_entries = []

for gene in list_of_genes:
    
    gene_label = f"*{gene.upper()}" if gene == "ha" else gene.upper()
    align_path = f"data/alignments/h3nx_{gene}.fasta"
    
    if gene == "ha":
        backbone_path = "data/ha/output.nwk"
        descriptor_entries.append([gene_label, align_path, backbone_path])
            
    else:
        div_path = f"results/trees_rooted/h3nx_{gene}_rooted/rerooted.newick"
        descriptor_entries.append([gene_label, align_path, div_path])

# descriptor csv
with open(f"{args.descriptor}", 'w') as descriptor_file:
	for row in descriptor_entries:
		descriptor_file.write(','.join(row) + '\n')
    
