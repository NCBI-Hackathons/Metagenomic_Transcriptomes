# -*- coding: utf-8 -*-
""" Utility functions for MetaTate package


"""
__author__ = "Matthew Bendall"

def merge_interval_list(ivs, dist=0):
    """ Merge intervals

    Args:
        ivs (list): List of intervals. Each interval is represented by a tuple of
            integers (start, end) where end > start.
        dist (int): Distance between intervals to be merged. Setting dist=1 will merge
            adjacent intervals

    Returns:
        list: Merged list of intervals

    Examples:
        >>> merge_intervals([])
        []
        >>> merge_intervals([(1,10)])
        [(1, 10)]
        >>> merge_intervals([(4, 9), (10, 14), (1, 3)])
        [(1, 3), (4, 9), (10, 14)]
        >>> merge_intervals([(4, 9), (10, 14), (1, 3)], dist=1)
        [(1, 14)]
    """
    if len(ivs)<= 1: return ivs
    ivs.sort(key=lambda x:x[0])
    ret = [ivs[0]]
    for iv in ivs[1:]:
        if iv[0] - ret[-1][1] > dist:
            ret.append(iv)
        else:
           ret[-1] = (ret[-1][0], max(iv[1],ret[-1][1]))
    return ret

def read_fasta(fh):
    ''' Reads multi-fasta data from file-like object
            fh  - file-like object with fasta data
        Generator yields tuples of strings (sequence_name, sequence)
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
    ''' Splits sequence onto multiple lines
            lwidth  - line width
    '''
    return '\n'.join(seq[i:i+lwidth] for i in range(0,len(seq), lwidth))
