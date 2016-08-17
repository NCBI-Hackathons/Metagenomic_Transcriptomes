#!/usr/bin/sh
## Script to generate a blast db for the COG protein sequences

cd databases/COG2014

gunzip -k prot2003-2014.fa.gz
makeblastdb -in prot2003-2014.fa -dbtype prot -out cog2003-2014

rm prot2003-2014.fa
