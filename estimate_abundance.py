#!/usr/bin/env python

'''
This will use bowtie2 and RSEM to estimate abundance for each individual COG00000_Trinity.fasta given 
a directory of Trinity output directories starting with COG
'''
import sys, os
 
def rsemSE(PATH, transcripts, reads, outprefix):
    cmd = 'perl '+PATH+' --transcripts '+transcripts+' --seqType fq --single '+reads
    cmd += ' --est_method RSEM --aln_method bowtie2 --trinity_mode --prep_reference --output_dir . --output_prefix '+outprefix
    os.system(cmd)
    
def rsemPE(PATH, transcripts, left_reads, right_reads, outprefix):    
    cmd = 'perl '+PATH+' --transcripts '+transcripts+' --seqType fq --left '+left_reads+' --right '+right_reads
    cmd += ' --est_method RSEM --aln_method bowtie2 --trinity_mode --prep_reference --output_dir . --output_prefix '+outprefix
    os.system(cmd)
        
if __name__ == "__main__":
    if len(sys.argv) == 4:
        DIR, PATH, reads = sys.argv[1:]
        for dir in os.listdir(DIR):
            if dir.endswith('.trinity') and os.path.isdir(dir):
                fasta = str(DIR)+'/'+dir+'/Trinity.fasta'
                outprefix = dir.split('.')[1]
                rsemSE(PATH, transcripts=fasta, reads=reads, outprefix=outprefix)
    elif len(sys.argv) == 5:
        DIR, PATH, left_reads, right_reads = sys.argv[1:]
        for dir in os.listdir(DIR):
            if dir.endswith('.trinity') and os.path.isdir(dir):
                fasta = str(DIR)+'/'+dir+'/Trinity.fasta'
                outprefix = dir.split('.')[1]
                rsemPE(PATH, transcripts=fasta, left_reads=left_reads, right_reads=right_reads, outprefix=outprefix)
    else:
        print "usage: python estimate_abundance.py <DIR with Trinity output directories> <PATH to Trinity 'align_and_estimate_abundance.pl'> <fastq file/s>"
        sys.exit()