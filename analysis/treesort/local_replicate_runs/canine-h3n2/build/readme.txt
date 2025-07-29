courtesy of jordan ort

in order to run treesort via snakemake you must do this:

git clone https://github.com/flu-crew/TreeSort.git
cd TreeSort
conda create -n nextstrain-treesort \
      --override-channels --strict-channel-priority \
      -c conda-forge -c bioconda --yes \
      augur auspice nextclade \
      snakemake git epiweeks \
      ncbi-datasets-cli csvtk seqkit tsv-utils \
      --file conda-requirements.txt
conda activate nextstrain-treesort
pip install .
nextstrain setup ambient
nextstrain build --ambient <YOUR_BUILD>
