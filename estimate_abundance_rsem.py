#!/usr/bin/env python

'''
This will use bowtie2 and RSEM to estimate abundance for each individual SRR0000.COG00000.assembled.fasta given 
a directory of Trinity output files
Requires Biopython, bowtie2 and rsem in your path
'''
import sys, os
import get_gene_trans_map

 
def rsemSE(transcripts, reads, outprefix):
    build = 'bowtie2-build '+transcripts+' '+transcripts+'.bowtie2'
    prep = 'rsem-prepare-reference  --transcript-to-gene-map '+transcripts+'.gene_trans_map '+transcripts+' '+transcripts+'.RSEM'
    bowtie = 'bowtie2 --no-mixed --no-discordant --gbar 1000 --end-to-end -k 200  -q -x '+transcripts+'.bowtie2 -U '+reads+' -p 4 | samtools view -F 4 -S -b -o '+outprefix+'.bowtie2.bam -'
    estimate = 'rsem-calculate-expression     -p 4 --fragment-length-mean 200 --fragment-length-sd 80   --no-bam-output --bam '+outprefix+'.bowtie2.bam '+transcripts+'.RSEM '+outprefix
    os.system(build)
    os.system(prep)
    os.system(bowtie)
    os.system(estimate)    
    os.system('rm *bowtie*')
    os.system('rm *RSEM*')
    
def rsemPE(transcripts, left_reads, right_reads, outprefix): 
    build = 'bowtie2-build '+transcripts+' '+transcripts+'.bowtie2'
    prep = 'rsem-prepare-reference  --transcript-to-gene-map '+transcripts+'.gene_trans_map '+transcripts+' '+transcripts+'.RSEM'
    bowtie = 'bowtie2 --no-mixed --no-discordant --gbar 1000 --end-to-end -k 200  -q -X 800 -x '+transcripts+'.bowtie2 -1 '+left_reads+' -2 '+right_reads+' -p 4 | samtools view -F 4 -S -b -o '+outprefix+'.bowtie2.bam -'
    estimate = 'rsem-calculate-expression  --paired-end   -p 4    --no-bam-output --bam '+outprefix+'.bowtie2.bam '+transcripts+'.RSEM '+outprefix
    os.system(build)
    os.system(prep)
    os.system(bowtie)
    os.system(estimate)    
    os.system('rm *bowtie*')
    os.system('rm *RSEM*')
    
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