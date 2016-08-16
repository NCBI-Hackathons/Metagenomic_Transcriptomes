#!/usr/bin/sh
## Generating simulated mock sequence data
## 	Mixture includes 40X E. coli, 20X Rhodo, and 10X Procholo
ART_PATH=$1

mkdir sim_seq
cd sim_seq

$1/art_illumina -sam -i ../data/ecoli/GCF_000012905.2_ASM1290v2_cds_from_genomic.fna -p -l 250 -ss MS -f 20 -m 750 -s 50 -o ecoli
$1/art_illumina -sam -i ../data/rhodo/GCF_000012905.2_ASM1290v2_cds_from_genomic.fna -p -l 250 -ss MS -f 40 -m 750 -s 50 -o rhodo
$1/art_illumina -sam -i ../data/prochloro/GCA_000011465.1_ASM1146v1_cds_from_genomic.fna -p -l 250 -ss MS -f 10 -m 750 -s 50 -o prochloro

# concatenating seqs
cat ecoli1.fq prochloro1.fq rhodo1.fq > sim_mock1.fq
cat ecoli2.fq prochloro2.fq rhodo2.fq > sim_mock2.fq
