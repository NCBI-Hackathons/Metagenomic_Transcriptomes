#!/usr/bin/env python

import sys
import subprocess
 
def blast(program, database, output, input=sys.stdin, evalue='0.01', nthreads='16'):
    cmd = '{prog} -query {input} -db {db} -max_target_seqs 1 -evalue {evalue} -outfmt 6 -out blastout -num_threads {nthreads}'.format(
            prog=program, input=input, out=output, db=database, evalue=evalue, nthreads=nthreads )
    blast_process = subprocess.Popen( cmd.split(), stdout=subprocess.PIPE )
    return_code = blast_process.wait()
    if return_code != 0:
        raise IOError("Blast Error!")
    return blast_process.stdout

if __name__ == "__main__":
    try:
        program, database, output, query_sequence_file, evalue, nthreads = sys.argv[1:]
    except ValueError:
        print 'Usage: blast_wrapper.py <program> <database> <output> <query_sequence_file> <evalue=0.01> <nthreads=16>'
        exit()
 
    result_handle = blast(program, database=database, output=output, input=query_sequence_file, evalue=evalue, nthreads=nthreads)
    print result_handle.read()
