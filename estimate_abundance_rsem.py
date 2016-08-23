#!/usr/bin/env python

'''
This will use bowtie2 and RSEM to estimate abundance for each individual SRR0000.COG00000.assembled.fasta given 
a directory of Trinity output files
Requires Biopython, bowtie2 and rsem in your path
'''
import sys, os
import get_gene_trans_map

 
def rsemSE(transcripts, reads, outprefix):
    prep = 'rsem-prepare-reference --bowtie2 --transcript-to-gene-map '+transcripts+'.gene_trans_map '+transcripts+' '+outprefix
    estimate = 'rsem-calculate-expression --bowtie2 '+reads+' '+outprefix+' '+outprefix
    os.system(prep)
    os.system(estimate)
    
def rsemPE(PATH, transcripts, left_reads, right_reads, outprefix):    
    prep = 'rsem-prepare-reference --bowtie2 --transcript-to-gene-map '+transcripts+'.gene_trans_map '+transcripts+' '+outprefix
    estimate = 'rsem-calculate-expression --bowtie2 --paired_end '+left_reads+' '+right_reads+' '+outprefix+' '+outprefix
    os.system(prep)
    os.system(estimate)    
    
if __name__ == "__main__":
    if len(sys.argv) == 3:
        DIR, reads = sys.argv[1:]
        for assembly in os.listdir(DIR):
            if assembly.endswith('.assembled.fasta'):
 #               fasta = str(DIR)+'/'+dir+'/Trinity.fasta'
                outprefix = assembly.split('.')[0]+'.'+assembly.split('.')[1]
                get_gene_trans_map.gene_trans(assembly, assembly+'.gene_trans_map')
                rsemSE(transcripts=DIR+'/'+assembly, reads=reads, outprefix=outprefix)
    elif len(sys.argv) == 4:
        DIR, left_reads, right_reads = sys.argv[1:]
        for assembly in os.listdir(DIR):
            if assembly.endswith('.assembled.fasta'):
#                fasta = str(DIR)+'/'+dir+'/Trinity.fasta'
                outprefix = assembly.split('.')[0]+'.'+assembly.split('.')[1]
                get_gene_trans_map.gene_trans(assembly, assembly+'.gene_trans_map')
                rsemPE(transcripts=DIR+'/'+assembly, left_reads=left_reads, right_reads=right_reads, outprefix=outprefix)
    else:
        print "usage: python estimate_abundance_rsem.py <DIR with Trinity output directories> <fastq file/s>"