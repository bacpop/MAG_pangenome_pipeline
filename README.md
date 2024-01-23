# MAG_pangenome_pipeline

A pipeline written in Snakemake to automatically generate pangenomes from metagenome assembled genomes (MAGs). 

## Dependencies: 

* Snakemake
* mmseqs2
* Bakta
* Biopython
* CheckM
* Pandas
* Rust

**NOTE:** Conda is used to call different environments and dependencies (see Snakemake file).

## To install:

Install the required packages using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)/[mamba](https://github.com/mamba-org/mamba):

```
git clone https://github.com/bacpop/MAG_pangenome_pipeline.git
mamba env create -f environment.yml
mamba activate celebrimbor
```

Download the required bakta database file:

```
bakta_db download --output /path/to/database
```

You can also use the light bakta database if using a suitable version of bakta:

```
bakta_db download --output /path/to/database --type light
```

Install [cgt](https://github.com/bacpop/cgt)

```
git clone https://github.com/bacpop/cgt.git
cd cgt
cargo build --release
```

## Quick start: 

Update `config.yaml` to specify workflow and directory paths. 
- `core`: gene frequency cutoff for core gene, anything above this frequency is annotated as a core gene.
- `output_dir`: path to output directory. Does not need to exist prior to running.
- `genome_fasta`: path to directory containing fasta files (must have either `.fa` or `.fasta` extension).
- `bakta_db`: path to bakta db downloaded above.
- `cgt_exe`: path to cgt executable. Relative path will be `cgt/target/release/cgt_rust`.
- `cgt_breaks`: frequency for rare/core gene cutoff, e.g. `0.1,0.9`, meaning genes predicted at `<0.1` frequency will be `rare`, `0.1<=x<0.9` will be `middle` and `>=0.9` will be `core`.
- `cgt_error`: sets false assignment rate of gene to particular frequency compartment.

Run snakemake (must be in same directory as `Snakemake` file):

```
snakemake -n <cores>
```

## Overview of workflow

This workflow annotates genes in metagenome-assembled genomes (MAGs) and using a probabilistic model to assign each gene to a gene frequency compartment based on their respective frequencies and genome completeness.

1. Predict genes in all FASTA files in given directory using [bakta](https://github.com/oschwengers/bakta)
1. Cluster genes using [mmseqs2](https://github.com/soedinglab/MMseqs2) and generate a gene presence/absence matrix
1. Generate a pangenome summary of observed gene frequencies
1. Calculate genome completeness using [CheckM](https://github.com/Ecogenomics/CheckM)
1. Probabistically assign each gene family as `core|middle|rare` using [cgt](https://github.com/bacpop/cgt)




