import re
import multiprocessing.pool as mp
import pandas

from collections import namedtuple


def chunks(l, n):
    """
    chunks a list into n-sized pieces
    """
    for i in range(0, len(l), n):
        # Create an index range for l of n items:
        yield l[i:i+n]


def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in seq[::-1]])


def parse_MD_tag(sequence, md_tag):
    tag = md_tag.split(':')[2]
    tag_chain = re.findall("([0-9]+)([AaCcGgTt]?)", tag)
    
    string = ''
    
    for link in tag_chain:
        dist, base = link
        
        string += '.' * int(dist)
        string += base
    
    return string

def geneset_to_pandas(geneset):
    """
    Life is a bit easier when a geneset is a pandas DataFrame.
    """
    items = []

    for n in geneset.dtype.names:
        v = geneset[n]
        # convert bytes columns to unicode (which pandas then converts to object)
        if v.dtype.kind == 'S':
            v = v.astype('U')
        items.append((n, v))

    return pandas.DataFrame.from_dict(dict(items))

Region = namedtuple('Region', ['chromosome', 'start', 'end', 'sequence'])