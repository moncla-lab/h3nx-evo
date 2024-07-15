#written by jordan ort in the moncla lab
#this removes any discrepancies introduced in the refine step, where sequences are taken out
#if they do not pass the clock filter. called in pruner step 

from Bio import SeqIO
from Bio import Phylo
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--alignments', type=str, required=True, nargs='+', help='path to alignment fasta files for each segment')
parser.add_argument('--trees', type=str, required=True, nargs='+', help='path to tree newick files for each segment') #timetree
parser.add_argument('--pruned_alignments', type=str, required=True, nargs='+', help='path to output alignment fasta files for each segment')
args = parser.parse_args()

tree_strains = {}

for tree in args.trees:
    t = Phylo.read(tree, 'newick')
    if not tree_strains:
        tree_strains = {strain:1 for strain in [clade.name for clade in t.find_clades() if 'NODE' not in clade.name]}
    else:
        for strain in [clade.name for clade in t.find_clades() if 'NODE' not in clade.name]:
            if strain in tree_strains:
                tree_strains[strain] += 1
            else:
                tree_strains[strain] = 1

unpruned_strains = [k for k,v in tree_strains.items() if v == len(args.trees)]

for fasta, outfasta in zip(args.alignments, args.pruned_alignments):
    records = []
    for record in SeqIO.parse(fasta, 'fasta'):
        if record.id in unpruned_strains:
            records.append(record)
        else:
            print(f'Pruning {record.id} as it only appeared in {tree_strains.get(record.id,0)} tree{"" if tree_strains.get(record.id,0)==1 else "s"}')
    SeqIO.write(records, outfasta, 'fasta')