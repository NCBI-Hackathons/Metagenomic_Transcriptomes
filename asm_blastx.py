#!/usr/bin/env python

''' 
Reciprocal BLAST of COG assemblies to COG protein sequences using blastx.
'''
import sys
import os
from blast_wrapper import blast

DIR = sys.argv[1]
BLASTDB=sys.argv[2]

for asm_dir in os.listdir(DIR):
	if "trinity" in asm_dir:
		print asm_dir
		asm_file = DIR + "/" + asm_dir + "/Trinity.fasta"
		blast_file = DIR + "/" + asm_dir + "/blastx.out"
		if os.path.isfile(asm_file): 
			blast(program="blastx",database=BLASTDB,input=asm_file, output=blast_file)
		else:
			print "No assembly file: " + asm_file

