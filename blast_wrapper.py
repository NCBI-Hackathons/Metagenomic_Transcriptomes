#!/usr/bin/env python

import sys
import subprocess
 
def blast(program, database, input=sys.stdin, evalue='0.01'):
    cmd = '{prog} -query {input} -db {db} -max_target_seqs 1 -evalue {evalue} -outfmt 6 -out blastout'.format(
            prog=program, input=input, db=database, evalue=evalue )
    blast_process = subprocess.Popen( cmd.split(), stdout=subprocess.PIPE )
    return_code = blast_process.wait()
    if return_code != 0:
        raise IOError("Blast Error!")
    return blast_process.stdout

if __name__ == "__main__":
    try:
        program, database, query_sequence_file, evalue = sys.argv[1:]
    except ValueError:
        print 'Usage: blast_wrapper.py <program> <database> <query_sequence_file> <evalue>'
        exit()
 
    result_handle = blast(program, database=database, input=query_sequence_file, evalue=evalue)
    print result_handle.read()