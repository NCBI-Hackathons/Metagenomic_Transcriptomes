# Metagenomic_Transcriptomes
Python pipeline to identify and functionally annotate expressed genes in an environment

## Purpose:
The advent of high-throughput sequencing methods has delivered a wealth of information to the scientific community, while simultaneously posing many new questions and creating significant data processing bottlenecks. One field high potential field for data-intensive applications is metagenomics, the culture-independent study of microbial genetic information from an environmental sample. Here we present an open-source, user-friendly software package for metagenomic functional analysis which bypasses taxonomic classification. The package utilizes existing software tools SRA tBLASTn, Trinity sequence assembly, RSEM transcript quantification, and BLASTx to create a streamlined pipeline for functional annotation and validation which can be used via the command line. Important features _____________. The pipeline was tested and validated using datasets ______________. We envision this software aiding to researchers seeking a functional understanding of metagenomic samples.

## Dependencies:
Python 2.7

## Directions to acquire code:
Clone the repository to a directory on your computer

Install MetaTate

```bash
git clone https://github.com/NCBI-Hackathons/Metagenomic_Transcriptomes.git
export PATH=$PWD/Metagenomic_Transcriptomes/bin:$PATH
```

Install dependencies

### `tblastn_vdb` (SRA Toolkit)

```bash
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.7.0/sratoolkit.2.7.0-centos_linux64.tar.gz
tar xzf sratoolkit.2.7.0-centos_linux64.tar.gz
export PATH=$PWD/sratoolkit.2.7.0-centos_linux64/bin:$PATH
```

```
build_COG_db


## Example run command:


## Estimated runtime:


