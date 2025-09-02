# PCRtaggingProfiler

A computational pipeline for characterizing and quantifying the genomic integration outcomes following PCR tagging.

## Project structure
- [data/](#-data)         --> Raw data (fastq.gz)
- [inputs/](#-inputs)     --> Barcode information, metadata, reference genome (human ref genome -- hg38.ucsc.fa -- is not uploaded)
- [envs/](#-envs)         --> Conda and Julia environment specification
- [julia/](#-julia)       --> Julia scripts for demultiplexing and trimming
- [R/](#-R)               --> R scripts for statistics
- [python](#-python)      --> Python scripts for mapping and filtering
- [Snakefile](#-snakefile) --> Main Snakemake workflow file
- [config.yml](#-config)   --> Configuration file for the pipeline

