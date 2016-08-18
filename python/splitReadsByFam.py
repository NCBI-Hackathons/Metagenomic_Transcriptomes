#! /usr/bin/env python
import sys,os
from collections import defaultdict
from itertools import chain
from glob import glob

def tab_line_gen(infile):
    ''' Returns a generator for tab-delmited file '''
    return (l.strip('\n').split('\t') for l in infile if not l.startswith('#'))

def parse_hitfile(fh):
    ret = defaultdict(set)
    for l in tab_line_gen(fh):
        ret[l[1]].add(l[0])
    return ret

def merge_dict_sets(dicts):
    ret = defaultdict(set)
    newkeys = set(list(chain.from_iterable(d.keys() for d in dicts)))
    for k in newkeys:
        ret[k] = set.union(*[d[k] for d in dicts])
    return ret

def fastqGen(fh):
    ''' Generator function for FASTQ file '''
    ret = []    
    iter = (l.strip('\n') for l in fh)
    for l in iter:
        ret.append(l)
        if len(ret)==4:
           yield ret
           ret = []

def main(args):
    dlist = [parse_hitfile(open(fn,'rU')) for fn in args.hits]
    # hitfiles = ['hits.short_1.txt', 'hits.short_3.txt']
    # d1 = parse_hitfile(open(hitfiles[0],'rU'))
    # d2 = parse_hitfile(open(hitfiles[1],'rU'))
    # Dictionary mapping family to spots
    fam_spots = merge_dict_sets(dlist)
    # Dictionary mapping spot to families
    spot_fams = defaultdict(set)
    for fam,spots in fam_spots.iteritems():
        for spot in spots:
            spot_fams[spot].add(fam) 
    
    # Remove files with prefix (because we are using append)
    for f in glob('%s.*.fq' % args.prefix):
        os.remove(f)
    
    for seq_fn in args.reads:
        iter = fastqGen(open(seq_fn,'rU'))
        for s in iter:
            spotid = s[0].split()[0].strip('@')
            for fam in spot_fams[spotid]:
                outfn = '%s.%s.fq' % (args.prefix, fam.split('|')[-1])
                with open(outfn, 'a') as outh:
                    print >>outh, '\n'.join(s)

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Split reads by family')
    parser.add_argument('--prefix', default='seqs/tmp',
                        help='''Prefix for output files''')        
    parser.add_argument('--hits', nargs='+',
                        help='''Output from BLAST''')    
    parser.add_argument('--reads', nargs='+',
                        help='''Output from BLAST''')    
    main(parser.parse_args())
