__author__ = 'Matthew L. Bendall'

def tab_line_gen(infile):
    ''' Returns a generator for tab-delmited file '''
    return (l.strip('\n').split('\t') for l in infile if not l.startswith('#'))

    


def fastqGen(fh):
    ''' Generator function for FASTQ file '''
    ret = []    
    iter = (l.strip('\n') for l in fh)
    for l in iter:
        ret.append(l)
        if len(ret)==4:
           yield ret
           ret = []

        yield list(itertools.islice(iter,4))

iter = (l.strip('\n') for l in open('short_1.fastq','rU'))

fq = fastqGen(open('short_1.fastq','rU'))
for rec in iter:
    pass
    
import itertools
itertools.islice(
        