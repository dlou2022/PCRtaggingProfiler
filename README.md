# PCRtaggingProfiler

A computational pipeline for characterizing and quantifying the genomic integration outcomes following PCR tagging.

## Project structure
data/           # Raw data (fastq.gz)
inputs/         # Barcode information, metadata, reference genome (human ref genome -- hg38.ucsc.fa -- is not uploaded)
envs/           # Conda and Julia environment specification
julia/          # Julia scripts for demultiplexing and trimming
R/              # R scripts
python/         # Python scripts
Snakefile       # Main Snakemake workflow file
config.yml      # Configuration file for the pipeline
README.md
