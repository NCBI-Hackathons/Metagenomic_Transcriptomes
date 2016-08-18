#! /usr/bin/env python

def main(args):
    print >>sys.stderr, args
    ''' Trinity --left --right '''

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Assemble')
    parser.add_argument('--assembler', default='trinity', choices=['trinity'],
                        help='''Assembler to use for family-wise assemblies.''')
    parser.add_argument('--SS_lib_type', 
                        help='''Strand-specific RNA-Seq read orientation for Trinity.
                                if paired: RF or FR, if single: F or R.''')
    parser.add_argument('--min_contig_length', type=int, default=200,
                        help=''' ''')
    parser.add_argument('--trimmomatic')
    parser.add_argument('--normalize_reads')
    parser.add_argument('--output',
                        help='''name of directory for output (will be created if it
                                doesn't already exist)''')
    parser.add_argument('hits', type=argparse.FileType('rU'),
                        help='''Hits file output by BLAST''')                       
    parser.add_argument('flat', type=argparse.FileType('rU'),
                        help='''Flat file to lookup references''')
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help="Output GTF")
    main(parser.parse_args())
