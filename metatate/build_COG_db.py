#! /usr/bin/env python
# -*- coding: utf-8 -*-
''' Build the COG database '''

__author__ = 'Matthew L. Bendall'


"""
#-----------------------------------------------------------
2.4.    cog2003-2014.csv

Contains list of orthology domains. Comma-delimited, format:

<domain-id>, <genome-name>, <protein-id>,<protein-length>,
<domain-start>, <domain-end>, <COG-id>, <membership-class>,

* Example:

333894695,Alteromonas_SN2_uid67349,333894695,427,1,427,COG0001,0,

* Comments:

In this version the fields <domain-id> and <protein-id> are identical
and both normally refer to GenBank GIs. Thus neither <domain-id> nor
<protein-id> are necessarily unique in this file (this happens when a
protein consists of more than one orthology domains, e.g. 48478501).

The <membership-class> field indicates the nature of the match
between the sequence and the COG consensus:

0 - the domain matches the COG consensus;

1 - the domain is significantly longer than the COG consensus;

2 - the domain is significantly shorter than the COG consensus;

3 - partial match between the domain and the COG consensus.

#-----------------------------------------------------------
"""
import sys
import re
import gzip
import urllib2
import os

from cStringIO import StringIO
from collections import defaultdict,Counter

from utils import merge_interval_list, read_fasta, fmt_seq

def main(args):
    # Mapping COG id to proteins
    cog_prot = defaultdict(lambda: defaultdict(list))
    # Mapping protein to COG ids
    prot_cog = defaultdict(lambda: defaultdict(list))

    if not os.path.exists(args.dest):
        print >>sys.stderr, 'Creating directory %s' % args.dest
        os.makedirs(args.dest)
    if not os.path.isdir(args.dest): sys.exit('Destination is not a directory')          
    
    if args.cog_csv:
        cog_csv = (l.strip('\n').split(',') for l in args.cog_csv)
    else:
        if not os.path.exists('%s/cog2003-2014.csv' % args.dest):
            print >>sys.stderr, 'Downloading COG csv from NCBI'
            response = urllib2.urlopen('ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cog2003-2014.csv')
            with open('%s/cog2003-2014.csv' % args.dest, 'w') as outh:
                outh.write(response.read())
        cog_csv = (l.strip('\n').split(',') for l in open('%s/cog2003-2014.csv' % args.dest, 'r'))
    
    for l in cog_csv:
        protid = l[2]
        cogid = l[6]
        spos = int(l[4])
        epos = int(l[5])
        cog_prot[cogid][protid] = merge_interval_list(cog_prot[cogid][protid] + [(spos,epos)], dist=1)
        prot_cog[protid][cogid] = merge_interval_list(prot_cog[protid][cogid] + [(spos,epos)], dist=1)
    
    print >>sys.stderr, "Found {:,} COGs".format(len(cog_prot))
    print >>sys.stderr, "Found {:,} protein ids".format(len(prot_cog))
    c = Counter(len(v) for k,v in prot_cog.iteritems())
    print >>sys.stderr, "Proteins belonging to multiple COGs:"
    print >>sys.stderr, 'Num COGs  | Count'
    for k in range(10):
        print >>sys.stderr, '%s%d' % (str(k+1).ljust(12),c[k+1])
    print >>sys.stderr,'%s%d' % ('11+'.ljust(12), sum(v for k,v in c.iteritems() if k>10))

    if args.fasta:
        seqiter = ((seqname,seq) for seqname,seq in read_fasta(gzip.GzipFile(fileobj=args.fasta)))
    else:
        if not os.path.exists('%s/prot2003-2014.fa.gz' % args.dest):
            print >>sys.stderr, 'Downloading protein sequences from NCBI'
            response = urllib2.urlopen('ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/prot2003-2014.fa.gz')
            with open('%s/prot2003-2014.fa.gz' % args.dest, 'wb') as outh:
                outh.write(response.read())
        seqiter = ((seqname,seq) for seqname,seq in read_fasta(gzip.open('%s/prot2003-2014.fa.gz' % args.dest, 'rb')))
    
    # Clear out destination
    purge_msg = True
    for cogid in cog_prot.keys():
        if os.path.exists('%s/%s/%s.faa' % (args.dest, cogid[:5], cogid)):
            if purge_msg:
                print >>sys.stderr, 'Purging existing sequence files'
                purge_msg = False
            os.remove('%s/%s/%s.faa' % (args.dest, cogid[:5], cogid))
    
    # Create the COG files
    numseqs = 0    
    filenames = {}
    seqcounts = Counter()
    for seqname,seq in seqiter:
        try:
            m = re.search('gi\|(\d+)\|ref', seqname)
            protid = m.group(1)
        except (AttributeError, ValueError):
            print >>sys.stderr, 'Error parsing sequence name: "%s"' % seqname
            continue
        numseqs += 1
        if not numseqs % 100000: print >>sys.stderr, 'Processed %d proteins...' % numseqs
        for cogid,ivs in prot_cog[protid].iteritems():
            for istart,iend in ivs:
                newseqname = 'cog|%s|%s (%d-%d)' % (cogid, seqname, istart, iend)
                newseq = seq[istart-1:iend]
                seqcounts[cogid] += 1
                if not os.path.exists('%s/%s' % (args.dest, cogid[:5])):
                    os.makedirs('%s/%s' % (args.dest, cogid[:5]))
                with open('%s/%s/%s.faa' % (args.dest, cogid[:5], cogid), 'a') as outh:
                    print >>outh, '>%s\n%s' % (newseqname,fmt_seq(newseq))
                filenames[cogid] = '%s/%s.faa' % (cogid[:5],cogid)
    
    print >>sys.stderr, "Processed %d total proteins" % numseqs
    with open('%s/cogfiles.txt' % args.dest, 'w') as outh:
        for k in sorted(filenames.keys()):
            print >>outh, '%s\t%s\t%d' % (k, filenames[k], seqcounts[k])


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Build the COG database.')
    parser.add_argument('--cog_csv', type=argparse.FileType('r'),
                        help='''Path to COG csv file (cog2003-2014.csv). Default is to download from NCBI.''')
    parser.add_argument('--fasta', type=argparse.FileType('rb'),
                        help='''Path to protein fasta file (prot2003-2014.fa.gz). Default is to download from NCBI.''')
    parser.add_argument('--dest', default='.', help="Destination directory.")
    main(parser.parse_args())
    
    