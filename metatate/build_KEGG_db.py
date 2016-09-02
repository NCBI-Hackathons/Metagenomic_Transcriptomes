#! /usr/bin/env python

import sys, re, time
import urllib2

# This is the number of records to fetch from NCBI (using efetch)
FETCHSIZE = 100

def attempt_url_open(url, max_attempts=5, quiet=False):
    ''' Attempt to open URL multiple times.
        This is useful in cases when urllib2.urlopen may fail randomly. For example, NCBI
        eutils sometimes fail when submitting an efetch query if it occurs too soon after
        the esearch query.
            max_attempts    - maximum number of attempts
            quiet           - if False, print information about failed attempts
        Returns file-like object if successful and None if URL could not be opened
    '''
    attempt = 1
    while attempt <= max_attempts:
        try:
            response = urllib2.urlopen(url)
            return response
        except urllib2.HTTPError as err:
            if not quiet:
                print >>sys.stderr, '        >>> Failed attempt %d' % attempt
                print >>sys.stderr, '        >>> HTTP Error: %s' % err.code
                print >>sys.stderr, '        >>> URL: %s' % url
            time.sleep(5)
    print >>sys.stderr, '        >>> Max attempts exceeded'
    print >>sys.stderr, '        >>> Failed to open URL: %s' % url    
    return None



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

def protein_accn_to_protein_seq(accn_list):
    ''' Retrieve protein sequences for accession or accession list
            accn_list   - list of accession numbers
            chunksize   - number of accessions to process for each request
        Returns generator that yields (sequence_name, sequence)
    '''
    if type(accn_list) is str:
        accn_list = [ accn_list ]
    if len(accn_list) > FETCHSIZE:
        print >>sys.stderr, "WARNING: Requesting %d accessions (Recommended size is <=%d)." % (len(accn_list), FETCHSIZE)
            
    search_args = {'term': '+OR+'.join('%s[accn]' % accn for accn in accn_list),
                   'usehistory': 'y',
                  }
    search_url = eutil_url('esearch', 'protein', search_args)
    
    fetch_args = None
    ncbi_error = 5
    while ncbi_error > 0:    
        response = attempt_url_open(search_url)
        if response is None:
            return
        data = response.read()
        m = re.search('<ERROR>(.+)</ERROR>', data)
        if m:
            ncbi_error = ncbi_error - 1
            continue
        else:
            try:
                fetch_args = {'rettype':   'fasta',
                              'retmode':   'text',
                              'WebEnv':    re.search('<WebEnv>(\S+)</WebEnv>', data).group(1),
                              'query_key': re.search('<QueryKey>(\d+)</QueryKey>', data).group(1),
                              }
                ncbi_error = 0
            except AttributeError:
                ncbi_error = ncbi_error - 1
    
    # Fetch sequences
    assert fetch_args is not None
    fetch_url = eutil_url('efetch', 'protein', fetch_args)
    response = attempt_url_open(fetch_url)
    if response is None:
        return
    else:
        return read_fasta(response)

""" Functions for interacting with KEGG API """
def all_KOs():
    ''' Retrieve all KEGG Orthology functions from KEGG
        Returns a generator expression that yields (ko_id, description)
    '''
    response = attempt_url_open('http://rest.kegg.jp/list/ko')
    if response is None:
        return
    return (tuple(l.strip('\n').split('\t')) for l in response if l.strip())

def kgenes_for_KO(ko_id):
    ''' Retrieve all KEGG genes for ortholog group
            ko_id   - KEGG ortholog ID, i.e. ko:K00162
        Returns generator that yields KEGG gene IDs
    '''
    ko_id = ko_id if ko_id.startswith('ko:') else 'ko:%s' % ko_id
    response = attempt_url_open('http://rest.kegg.jp/link/genes/%s' % ko_id)
    if response is None:
        return
    for l in response:
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
        response = attempt_url_open(url)
        if response is None:
            continue
        for l in response:
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
            ko_desc = list(all_KOs())
            if not ko_desc: sys.exit('>>> Failed to download KO ids')
            for ko_id, description in ko_desc:
                ko_ids.append(ko_id)
                print >>outh, '%s\t%s' % (ko_id, description)
        print >>sys.stderr, 'Downloaded %d KO ids, saved to %s' % (len(ko_ids), ko_txt)
    
    # Simplify KO ids
    tmp = []
    for ko_id in ko_ids:
        m = re.match('(ko:)?(K\d{5})', ko_id)
        if m: tmp.append(m.group(2))
        else: print >>sys.stderr, 'Improperly formatted ID: %s' % ko_id
    ko_ids = tmp
    
    # User supplied KO ids
    if args.ko_ids is not None:
        tmp = []
        for k in args.ko_ids.split(','):
            m = re.match('(ko:)?(K\d{5})', k)
            # Check that id is formatted correctly and exists in KEGG
            if m and m.group(2) in ko_ids:
                tmp.append(m.group(2))
            else: print >>sys.stderr, 'Unknown ID: %s' % k
        ko_ids = tmp
    if not ko_ids:
        sys.exit('No KO ids found')

    for ko_id in ko_ids:
        print >>sys.stderr, '> Processing %s...' % ko_id
        if args.use_subdirs:
            kodir = '%s/%s/%s' % (args.dest, ko_id[:4], ko_id)
            if not os.path.exists(kodir): os.makedirs(args.dest)
            gene_txt = '%s/%s.genes.txt' % (kodir, ko_id)
        else:
            gene_txt = '%s/%s.genes.txt' % (args.dest, ko_id)
        if os.path.exists(gene_txt):
            accn_list = [l.strip('\n').split('\t')[1] for l in open(gene_txt, 'rU')]
            print >>sys.stderr, '    Found %d accessions in %s' % (len(accn_list), gene_txt)
        else:
            print >>sys.stderr, '    Downloading gene list for %s' % ko_id
            kgene_ids = list(kgenes_for_KO(ko_id))
            if not kgene_ids:
                print >>sys.stderr, '    >>> Failed to download gene list for %s' % ko_id
                continue
            with open(gene_txt, 'w') as outh:
                accn_list = []
                for kgid,accn in ncbi_protein_accn_for_kgenes(kgene_ids):
                    print >>outh, '%s\t%s' % (kgid,accn)
                    accn_list.append(accn)
            print >>sys.stderr, '    Downloaded %d accessions, saved to %s' % (len(accn_list), gene_txt)
        
        # Download sequences
        if args.use_subdirs:
            kodir = '%s/%s/%s' % (args.dest, ko_id[:4], ko_id)
            if not os.path.exists(kodir): os.makedirs(args.dest)
            gene_faa = '%s/%s.faa' % (kodir, ko_id)
        else:
            gene_faa = '%s/%s.faa' % (args.dest, ko_id.split(':')[1])
        seqcount = 0
        with open(gene_faa, 'w') as outh:
            for i in range(0, len(accn_list), FETCHSIZE):
                accn_subset = accn_list[i:i+FETCHSIZE]
                iter = protein_accn_to_protein_seq(accn_subset)
                if iter is None:
                    print >>sys.stderr, '    >>> Sequence download failed for %d sequences.\n%s' % (len(accn_subset), ','.join(accn_subset))
                    continue
                for n,s in iter:
                    seqcount += 1
                    print >>outh, '>%s\n%s' % (n, fmt_seq(s))
        
        print >>sys.stderr, '    Wrote %d sequences to %s' % (seqcount, gene_faa)

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Build KEGG database')
    parser.add_argument('--use_subdirs', action='store_true',
                        help='''Create subdirectories within destination directory. The
                                full KEGG database has over 20K KOs. This many files does
                                not play nicely with some filesystems, so you may want to
                                break it up.'''
                        )    
    parser.add_argument('--ko_ids',
                        help='''Comma seperated list of KO ids to download. Default is to 
                                download all.''')
    parser.add_argument('dest', default='.', nargs='?',
                        help='''Destination directory''')    
    main(parser.parse_args())



