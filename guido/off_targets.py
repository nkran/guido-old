import os
import subprocess
import tempfile
import pandas as pd
from io import StringIO
from tqdm import tqdm

import guido.log as log
from guido.helpers import chunks, rev_comp

logger = log.createCustomLogger('off-targets')
ROOT_PATH = os.path.abspath(os.path.dirname(__file__))


def run_bowtie(cut_sites, max_offtargets, threads):
    genome_index_path = os.path.join(ROOT_PATH, 'data', 'references', 'AgamP4')
    mismatches = 3

    # create temporary file for bowtie input
    temp = tempfile.NamedTemporaryFile(mode='w+t', prefix="guido_")
    temp.write('\n'.join(['>seq|{}|{}|{}|{}\n{}'.format(i, cut['guide_loc'][0], cut['guide_loc'][1], cut['guide'], cut['guide']) for i, cut in enumerate(cut_sites)]))

    # run bowtie alignment
    bowtie_command = 'bowtie -p {} -v {} --sam --sam-nohead -k {} {} -f {}'.format(threads, mismatches, max_offtargets, genome_index_path, temp.name)
    rproc = subprocess.Popen(bowtie_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    stdout, stderr = rproc.communicate()

    # output bowtie to pandas dataframe
    targets = pd.read_csv(StringIO(stdout), sep='\t', names=list(range(14)),
                          header=None, index_col=False, converters={3:int})
    
    # split sequence identification
    nsplit = targets[0].str.split('|', n = 4, expand = True)
    nsplit.columns = ['s', 'id', 'start', 'end', 'seq']

    # add it back to the target df
    targets = pd.concat([targets, nsplit.iloc[:,1:]], axis=1)

    # extract the number of mismatches in the off-target
    targets['mm'] = targets.apply(lambda x: x[13].split(':')[-1], axis=1)

    # count off-targets for a given guide
    mismatches_count = targets.groupby(['id'])['mm'].value_counts().unstack().fillna(0).reset_index()
    mismatches_count['id'] = mismatches_count['id'].astype(int)
    
    # add mismatch info to the guide dict
    for ix, x in mismatches_count.iterrows():
        i = int(x.id)
        cut_sites[i]['mm'] = x[['0', '1', '2', '3']].to_dict()
    
    return cut_sites, targets

