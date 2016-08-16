#!/usr/bin/env python

'''
usage: python rpstblastn_wrapper.py <transcripts fasta> <path to cog index> <evalue cutoff> <nthreads>
'''
import sys, os

def rpstblastn(transcripts,cog,evalue,nthreads):
    output = str(transcripts).split('/')[-1]+'.rpstblastnout'
    print output
    cmd = 'rpstblastn -query '+transcripts+' -db '+cog+' -num_threads '+nthreads+' -evalue '+evalue+' -out '+output+' -outfmt 6 -max_target_seqs 1'
    os.system(cmd)
    
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print "python rpstblastn_wrapper.py <transcripts fasta> <path to cog index> <evalue cutoff> <nthreads>"
        sys.exit(0)

    transcripts,cog,evalue,nthreads = sys.argv[1:]
    rpstblastn(transcripts,cog,evalue,nthreads)
