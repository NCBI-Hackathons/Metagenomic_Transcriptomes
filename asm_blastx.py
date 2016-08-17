
#!/usr/bin/env python

''' 
Reciprocal BLAST of COG assemblies to COG protein sequences using blastx.
'''
import sys
import os
from blast_wrapper import blast

DIR = sys.argv[1]
BLASTDB=sys.argv[2]

for asm_file in os.listdir(DIR):
	asm_file = DIR + "/" + asm_file 
	if "assembled.fasta" in asm_file:
		blast_file = asm_file.replace("assembled.fasta", "assembled.blastx.tsv")
		blast(program="blastx",database=BLASTDB,input=asm_file, output=blast_file)
