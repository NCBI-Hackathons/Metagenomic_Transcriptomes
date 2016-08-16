#!/usr/bin/env python

'''
usage: python rpstblastn_wrapper.py <SRA fasta> <path to cog index> <evalue cutoff> <nthreads>
'''
import sys, os

def rpstblastn(SRA,cog,evalue,nthreads):
    output = str(SRA).split('/')[-1]+'.rpstblastnout'
    print output
    cmd = 'rpstblastn -query '+SRA+' -db '+cog+' -num_threads '+nthreads+' -evalue '+evalue+' -out '+output+' -outfmt 6 -max_target_seqs 1'
    os.system(cmd)
    
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print "python rpstblastn_wrapper.py <SRA fasta> <path to cog index> <evalue cutoff> <nthreads>"
        sys.exit(0)

    SRA,cog,evalue,nthreads = sys.argv[1:]
    rpstblastn(SRA,cog,evalue,nthreads)
