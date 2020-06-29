<p align="center">
<img src="Metasequencing.png" alt="metasequencing logo">
</p>
<br>

#### A pipeline for taxonomic/functional annotation and quantification of metagenomic sequencing data.

**README is a work in progress.**

MetaSequencing takes a partial assembly approach to processing metagenomic data, assembling short-reads to contigs that are used for taxonomic annotation and de novo open-reading frame detection and functional annotation.
This provides a trade off between short-read only methods and full metagenomic genome assembly, conceptually increasing taxonomic assignment accuracy by using longer contigs whilst negating the need to deal with binning and associated issues as when trying to reconstruct full genomes. This also enables quantification of function at the per taxon level.

The pipeline can be installed into isolated Conda environments using the attached install script and full metagenomic processing can be run on FASTA/FASTQ files using a single command, generating summary reports alongside processed data tables. MetaSequencing is built using the [SnakeMake](https://snakemake.readthedocs.io/en/stable/) workflow management system. A previous version of the pipeline based on the [CGAT](https://github.com/cgat-developers/cgat-core) workflow is available [here](https://github.com/microbialman/CGATMetaSequencing) but lacks some improvements made in this version and is no longer actively developed.

*Note: The complete pipeline has been tested in-house but no formal benchmarking is currently available. However, each of the methods used have been published and are widely used and hopefully some formal testing will be done soon!*

## Pipeline overview

*Details to follow.*

## Install

MetaSequencing and its dependencies are best installed and run within a self-contained Conda environment.
Following the instructions below will create the necessary environments and install all dependencies.
This assumes an existing working install of [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

```bash
#clone this repository into the folder where you would like the install
git clone https://github.com/microbialman/MetaSequencing.git
cd MetaSequencing
#install of the main metasequencing environment may take several minutes as it solves
sh install.sh
```

*Note: This will generate two conda environments (some tools require Python 2 in place of 3). If you are manually installing MetaSequencing to environments without the default names you will also need to edit the calls that load the Python 2 environment in the config .yaml files as needed*

## Running MetaSequencing

*Detailed documentation will be added soon.*
