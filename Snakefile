
"""Here, define your wildcards. To include more subtypes or gene segments, simply
add those to these lists, separated by commas"""
SUBTYPES = ["h3nx"]
SEGMENTS = ["ha","pb2","pb1","na","np","pa","ns","mp"]
#SEGMENTS = ["ha"]

"""This rule tells Snakemake that at the end of the pipeline, you should have
generated JSON files in the auspice folder for each subtype and segment."""
rule all:
    input:
        auspice_json = expand("auspice/nonhuman_influenza_{subtype}_{segment}.json", subtype=SUBTYPES, segment=SEGMENTS)

"""Specify all input files here. For this build, you'll start with input sequences
from the example_data folder, which contain metadata information in the
sequence header. Specify here files denoting specific strains to include or drop,
references sequences, and files for auspice visualization (colors)"""
rule files:
    params:
        input_sequences = "example_data/{subtype}_{segment}.fasta",
        dropped_strains = "config/dropped_strains_{subtype}.txt",
        include_strains = "config/include_strains_{subtype}.txt",
        reference = "config/references/{subtype}_{segment}.gb", #H3N8 from 1978
        auspice_config = "config/auspice_config_{subtype}.json",
        colors = "config/colors_{subtype}.tsv",
        lat_long = "config/lat_longs_{subtype}.tsv"

files = rules.files.params


"""The minimum length required for sequences. Sequences shorter than these will be
subsampled out of the build. Here, we're requiring all segments to be basically
complete. To include partial genomes, shorten these to your desired length"""

def group_by(w):
    gb = {'h3nx': 'year country host subtype'}
    return gb[w.subtype] 

def sequences_per_group(w):
    spg = {'h3nx':'30'}
    return spg[w.subtype] 

def min_length(w):
    len_dict = {"pb2": 2100, "pb1": 2100, "pa": 2000, "ha":1600, "np":1400, "na":1270, "mp":900, "ns":800}
    length = len_dict[w.segment]
    return(length)

"""Sequences with sample collection dates earlier than these will be subsampled out of the build"""
def min_date(w):
    date = {'h3nx':'1960'}
    return date[w.subtype]

def traits_columns(w):
    traits = {'h3nx':'host country region subtype order'}
    return traits[w.subtype]

"""In this section of the Snakefile, rules are specified for each step of the pipeline.
Each rule has inputs, outputs, parameters, and the specific text for the commands in
bash. Rules reference each other, so altering one rule may require changing another
if they depend on each other for inputs and outputs. Notes are included for
specific rules."""


"""The parse rule is used to separate out sequences and metadata into 2 distinct
files. This rule assumes an input fasta file that contains metadata information
in the header. By specifying the order of those fields in the `fasta_fields` line,
`augur parse` will separate those fields into labeled columns in the output metadata
file."""
rule parse:
    message: "Parsing fasta into sequences and metadata"
    input:
        sequences = files.input_sequences
    output:
        sequences = "results/sequences/sequences_{subtype}_{segment}.fasta",
        metadata = "results/metadata/metadata_{subtype}_{segment}.tsv"
    params:
        fasta_fields =  "strain accession subtype date host country region species broad order",
        prettify_fields = "country host species region order"
    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --fields {params.fasta_fields} \
            --prettify-fields {params.prettify_fields}
        """

"""This rule specifies how to subsample data for the build, which is highly
customizable based on your desired tree."""
rule filter:
    message:
        """
        Filtering to
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
          - excluding strains in {input.exclude}
          - samples with missing region and country metadata
          - excluding strains prior to {params.min_date}
        """
    input:
        sequences = rules.parse.output.sequences,
        metadata = rules.parse.output.metadata,
        exclude = files.dropped_strains,
        include = files.include_strains
    output:
        sequences = "results/filtered/filtered_{subtype}_{segment}.fasta"
    params:
        group_by = group_by, 
        sequences_per_group = sequences_per_group,
        min_date = min_date,
        min_length = min_length,
        exclude_where = "host=ferret"

    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --include {input.include} \
            --output {output.sequences} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --exclude-where {params.exclude_where} \
            --min-length {params.min_length} \
            --non-nucleotide
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = files.reference
    output:
        alignment = "results/alignments/aligned_{subtype}_{segment}.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --remove-reference \
            --nthreads 1
        """


rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/iqtree/raw_tree/tree-raw_{subtype}_{segment}.nwk"
    params:
        method = "iqtree"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method {params.method} \
            --nthreads 1
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = rules.parse.output.metadata
    output:
        tree = "results/iqtree/timed_tree/tree_{subtype}_{segment}.nwk",
        node_data = "results/iqtree/branch_lengths/branch-lengths_{subtype}_{segment}.json"
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_filter_iqd = 4
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/tree/nt_muts/nt-muts_{subtype}_{segment}.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}\
            --keep-ambiguous
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "results/tree/aa_muts/aa-muts_{subtype}_{segment}.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data}
        """

rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata
    output:
        node_data = "results/tree/traits/traits_{subtype}_{segment}.json",
    params:
        columns = traits_columns,
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

"""This rule exports the results of the pipeline into JSON format, which is required
for visualization in auspice. To make changes to the categories of metadata
that are colored, or how the data is visualized, alter the auspice_config files"""
rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata,
        node_data = [rules.refine.output.node_data,rules.traits.output.node_data,rules.ancestral.output.node_data,rules.translate.output.node_data],
        auspice_config = files.auspice_config,
	colors = files.colors,
	lat_long = files.lat_long
    output:
        auspice_json = "auspice/nonhuman_influenza_{subtype}_{segment}.json"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data}\
            --auspice-config {input.auspice_config} \
            --include-root-sequence \
			--colors {input.colors} \
			--lat-longs {input.lat_long} \
            --output {output.auspice_json}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
