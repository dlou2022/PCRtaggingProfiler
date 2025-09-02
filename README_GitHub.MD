# PCRtaggingProfiler

A computational pipeline for characterizing and quantifying the genomic integration outcomes following PCR tagging.

![License](https://img.shields.io/badge/license-MIT-blue.svg) ![Python Version](https://img.shields.io/badge/python-3.8%2B-brightgreen) ![Status](https://img.shields.io/badge/status-active-success)

---

## 📖 Project Structure
.
├── data/           # Raw data (fastq.gz)
├── inputs/         # Barcode information, metadata, reference genome (human ref genome -- hg38.ucsc.fa -- is not uploaded)
├── envs/           # Conda and Julia environment specification
├── julia/          # Julia scripts for demultiplexing and trimming
├── R/              # R scripts
├── python/         # Python scripts
│
├── Snakefile       # Main Snakemake workflow file
├── config.yml      # Configuration file for the pipeline
└── README.md


