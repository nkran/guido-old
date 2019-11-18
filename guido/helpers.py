import re
import multiprocessing.pool as mp

def istarmap(self, func, iterable, chunksize=1):
    """
    starmap-version of imap
    """
    if self._state != mp.RUN:
        raise ValueError("Pool not running")

    if chunksize < 1:
        raise ValueError(
            "Chunksize must be 1+, not {0:n}".format(
                chunksize))

    task_batches = mp.Pool._get_tasks(func, iterable, chunksize)
    result = mp.IMapIterator(self._cache)
    self._taskqueue.put(
        (
            self._guarded_task_generation(result._job,
                                          mp.starmapstar,
                                          task_batches),
            result._set_length
        ))
    return (item for chunk in result for item in chunk)

mp.Pool.istarmap = istarmap


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