#! /usr/bin/env python

import urllib2
import re

""" Functions for interacting with NCBI eutils """
def eutil_url(app, db, args):
    ''' Returns a formatted eutils URL
            app     - eutil application (i.e. esearch, efetch)
            db      - database to query
            args    - additional args provided to eutils
    '''
    baseurl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
    argstr = '&'.join('%s=%s' % (k,v) for k,v in args.iteritems())
    return '%(baseurl)s/%(app)s.fcgi?db=%(db)s&%(argstr)s' % \
        {'baseurl': baseurl, 'app': app, 'db': db, 'argstr': argstr}

def read_fasta(fh):
    ''' Reads multi-fasta data from file-like object
            fh  - file-like object with fasta data
        Generator yields tuples: (sequence_name, sequence)
    '''
    curname, curseq  = (None, '')
    for l in fh:
        if l.startswith('>'):
            if curname is not None:
                yield (curname, curseq)
            curname = l.strip('\n').lstrip('>')
            curseq = ''
        else:
            curseq += l.strip('\n')
    yield (curname, curseq)

def fmt_seq(seq, lwidth=70):
    ''' Returns multi-line sequence 
            lwidth  - line width
    '''
    return '\n'.join(seq[i:i+lwidth] for i in range(0,len(seq), lwidth))

def protein_accn_to_protein_seq(accn_list, chunksize=100):
    ''' Retrieve protein sequences for accession or accession list
            accn_list   - list of accession numbers
            chunksize   - number of accessions to process for each request
        Generator yields tuples: (sequence_name, sequence)
    '''
    if type(accn_list) is str:
        accn_list = [ accn_list ]
    for i in range(0, len(accn_list), chunksize):
        search_args = {'term': '+OR+'.join('%s[accn]' % accn for accn in accn_list[i:i+chunksize]),
                       'usehistory': 'y',
                      }
        search_url = eutil_url('esearch', 'protein', search_args)
        data = urllib2.urlopen(search_url).read()
        # Parse out the history tokens
        fetch_args = {'WebEnv':   re.search('<WebEnv>(\S+)</WebEnv>', data).group(1),
                      'query_key': re.search('<QueryKey>(\d+)</QueryKey>', data).group(1),
                      'rettype': 'fasta',
                      'retmode': 'text',
                     }
        fetch_url = eutil_url('efetch', 'protein', fetch_args)
        for sname,seq in read_fasta(urllib2.urlopen(fetch_url)):
            yield sname.strip('>'), seq

""" Functions for interacting with KEGG API """
def all_KOs():
    ''' Retrieve all KEGG Orthology functions from KEGG'''
    for l in urllib2.urlopen('http://rest.kegg.jp/list/ko'):
        if l:
            yield l.strip('\n').split('\t')

def kgenes_for_KO(ko_id):
    ''' Retrieve all KEGG genes for ortholog group
            ko_id   - KEGG ortholog ID, i.e. ko:K00162
        Returns generator that yields KEGG gene IDs
    '''
    ko_id = ko_id if ko_id.startswith('ko:') else 'ko:%s' % ko_id
    url = 'http://rest.kegg.jp/link/genes/%s' % ko_id
    for l in urllib2.urlopen(url):
        m = re.match('%s\t(\S+)' % ko_id, l)
        if m: yield m.group(1)

def ncbi_protein_accn_for_kgenes(kgene_ids, chunksize=10):
    ''' Retrieve NCBI accessions for KEGG genes
            kgene_ids    - List of KEGG gene ids, i.e. dpo:Dpse_GA17214
        Returns generator that yields (kgene_id, ncbi_acc)
    '''
    if type(kgene_ids) is str:
        kgene_ids = [ kgene_ids ]
    for i in range(0,len(kgene_ids),chunksize):
        url = 'http://rest.kegg.jp/conv/ncbi-proteinid/%s' % '+'.join(kgene_ids[i:i+chunksize])
        for l in urllib2.urlopen(url):
            m = re.match('(\S+)\tncbi-proteinid:(\S+)', l)
            if m: yield m.groups()

def main(args):
    import sys, os
    if not os.path.exists(args.dest):
        os.makedirs(args.dest)
    else:
        if not os.path.isdir(args.dest): sys.exit("Destination is not directory")
    
    assert os.path.exists(args.dest) and os.path.isdir(args.dest)
    
    # Identify KEGG Orthology IDs
    ko_txt = '%s/KEGG_Orthology.txt' % args.dest
    if os.path.exists(ko_txt):
        ko_ids = [l.strip('\n').split('\t')[0] for l in open(ko_txt, 'rU')]
        print >>sys.stderr, 'Found %d KO ids in %s' % (len(ko_ids), ko_txt)
    else:
        print >>sys.stderr, 'Downloading KO ids'    
        ko_ids = []
        with open(ko_txt, 'w') as outh:
            for ko_id, description in all_KOs():
                ko_ids.append(ko_id)
                print >>outh, '%s\t%s' % (ko_id, description)
        print >>sys.stderr, 'Downloaded %d KO ids, saved to %s' % (len(ko_ids), ko_txt)
    
    # User supplied KO ids
    if args.ko_ids is not None:
        tmp = [k if k.startswith('ko:') else 'ko:%s' % k for k in args.ko_ids.split(',')]
        for k in tmp:
            if k not in ko_ids: sys.exit('%s not found in KEGG' % k)
        ko_ids = tmp

    for ko_id in ko_ids:
        print >>sys.stderr, '> Processing %s...' % ko_id
        gene_txt = '%s/%s.genes.txt' % (args.dest, ko_id.split(':')[1])
        if os.path.exists(gene_txt):
            accn_list = [l.strip('\n').split('\t')[1] for l in open(gene_txt, 'rU')]
            print >>sys.stderr, '    Found %d accessions in %s' % (len(accn_list), gene_txt)
        else:
            print >>sys.stderr, '    Downloading gene list for %s' % ko_id        
            with open(gene_txt, 'w') as outh:
                kgene_ids = list(kgenes_for_KO(ko_id))
                accn_list = []
                for kgid,accn in ncbi_protein_accn_for_kgenes(kgene_ids):
                    print >>outh, '%s\t%s' % (kgid,accn)
                    accn_list.append(accn)
            print >>sys.stderr, '    Downloaded %d accessions, saved to %s' % (len(accn_list), gene_txt)
        
        # Download sequences
        gene_faa = '%s/%s.faa' % (args.dest, ko_id.split(':')[1])
        seqcount = 0
        with open(gene_faa, 'w') as outh:
            for n,s in protein_accn_to_protein_seq(accn_list):
                seqcount += 1
                print >>outh, '>%s\n%s' % (n, fmt_seq(s))
        print >>sys.stderr, '    Wrote %d sequences to %s' % (seqcount, gene_faa)

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Build KEGG database')
    parser.add_argument('--ko_ids',
                        help='''Comma seperated list of KO ids to download. Default is to 
                                download all.''')
    parser.add_argument('dest', default='.', nargs='?',
                        help='''Destination directory''')    
    main(parser.parse_args())



