
#!/usr/bin/env python

''' 
Check reciprical blast results to make sure protein blast hit matches assembly COG
'''
import sys
import os

IDTBL=sys.argv[1]
DIR=sys.argv[2]

with open(IDTBL, 'r') as tsv:
	cog_tbl = [line.strip().split('\t') for line in tsv]

## Need to check output then work on match checking
for blast_file in os.listdir(DIR):

	if "assembled.blastx.tsv" in blast_file:
		print blast_file
		blast_cog = blast_file.split(".")[1]
		print blast_cog
		blast_file = DIR + "/" + blast_file

		with open(blast_file,'r') as tsv:
			blast_out = [line.strip().split('\t') for line in tsv]

		for row in blast_out:
			print row	
