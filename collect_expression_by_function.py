#!/usr/bin/env python

'''
This will read all of the isoform.results files and add expression values to functions by COG
Take as input the COG # to function table and a directory with COG#.isoform.results files
'''
import sys, os
import pandas as pd
import numpy
import csv

#COG functional categories
A = []
B = []
C = []
D = []
E = []
F = []
G = []
H = []
I = []
J = []
K = []
L = []
M = []
N = []
O = []
P = []
Q = []
T = []
U = []
Y = []
Z = []
R = []
S = []

COG_map = sys.argv[1]
COG_function_dict = {}
COG = open(COG_map, 'rU')
for line in COG:
    num = line.split('\t')[0]
    fun = line.split('\t')[1]
    COG_function_dict[num] = fun

DIR = sys.argv[2]
for file in os.listdir(DIR):
    if file.endswith('.isoforms.results'):
        cog_num = file.split('.')[0]
        results = pd.read_csv(str(DIR)+'/'+file, sep='\t')
        count = sum(results['expected_count'])
        fun = COG_function_dict[cog_num]