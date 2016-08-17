#!/usr/bin/env python

'''
This will read all of the isoform.results files and add expression values to functions by COG
Take as input the COG # to function table and a directory with SRR#.COG#.isoform.results files
'''
import sys, os
import pandas as pd
import numpy
from collections import defaultdict

#COG functional categories
fun_dict = {
'A' : 0,
'B' : 0,
'C' : 0,
'D' : 0,
'E' : 0,
'F' : 0,
'G' : 0,
'H' : 0,
'I' : 0,
'J' : 0,
'K' : 0,
'L' : 0,
'M' : 0,
'N' : 0,
'O' : 0,
'P' : 0,
'Q' : 0,
'T' : 0,
'U' : 0,
'Y' : 0,
'Z' : 0,
'R' : 0,
'S' : 0
}
#create COG# function dictionary
COG_map = sys.argv[1]
COG_function_dict = {}
COG = open(COG_map, 'rU')
for line in COG:
    num = line.split('\t')[0]
    fun = line.split('\t')[1]
    COG_function_dict[num] = fun

#Get expected count for each COG assembly
DIR = sys.argv[2]
for file in os.listdir(DIR):
    if file.endswith('.isoforms.results'):
        cog_num = file.split('.')[1]
        results = pd.read_csv(str(DIR)+'/'+file, sep='\t')
        count = sum(results['expected_count'])
        fun = COG_function_dict[cog_num]
        fun_dict[fun] += count

print fun_dict