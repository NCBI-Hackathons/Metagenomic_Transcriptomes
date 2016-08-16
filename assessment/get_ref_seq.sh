#!/usr/bin/sh
## Downloading reference genome annotation data from RefSeq for three strain mock community
##	Mock community biosample - http://www.ncbi.nlm.nih.gov/bioproject/85559
##	Strains 
##	* Escherichia coli K-12 MG1655
##	* Rhodobacter sphaeroides 2.4.1
##	* Prochlorococcus marinus MED4

## Directory Structure
mkdir -p data/{ecoli,rhodo,prochloro}

## Downloading data
cd data/ecoli
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000012905.2_ASM1290v2/*

cd ../rhodo
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000012905.2_ASM1290v2/*

cd ../prochloro
wget -r ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000011465.1_ASM1146v1/*
