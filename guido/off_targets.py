import os
import subprocess
import tempfile
import pandas as pd
import numpy as np
from io import StringIO

import guido.log as log

logger = log.createCustomLogger('off-targets')


def run_bowtie(cut_sites, max_offtargets, genome_index_path, threads):
    mismatches = 3

    # create temporary file for bowtie input
    with tempfile.NamedTemporaryFile(mode='w+t', prefix="guido_") as temp:
        temp.write(
            '\n'.join(
                [
                    '>seq|{}|{}|{}|{}\n{}'.format(
                        i,
                        cut['guide_loc'][0],
                        cut['guide_loc'][1],
                        cut['guide'],
                        cut['guide'],
                    )
                    for i, cut in enumerate(cut_sites)
                ]
            )
        )
        temp.flush()
        # run bowtie alignment
        # bowtie_command = 'bowtie -p {} -v {} --sam --sam-nohead -k {} {} -f {}'.format(
        #     threads, mismatches, max_offtargets, genome_index_path, temp.name
        # )
        bowtie_command = 'bowtie -p {} -v {} --sam --sam-nohead -a {} -f {}'.format(
            threads, mismatches, genome_index_path, temp.name
        )
        rproc = subprocess.Popen(
            bowtie_command.split(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        stdout, stderr = rproc.communicate()

    # output bowtie to pandas dataframe
    targets = pd.read_csv(
        StringIO(stdout),
        sep='\t',
        names=list(range(14)),
        header=None,
        index_col=False,
        converters={3: int},
    )

    # split sequence identification
    nsplit = targets[0].str.split('|', n=4, expand=True)
    nsplit.columns = ['s', 'id', 'chrom', 'start', 'seq']

    # add it back to the target df
    targets = pd.concat([targets, nsplit.iloc[:, 1:]], axis=1)
    targets['id'] = targets['id'].astype(int)
    targets['start'] = targets['start'].astype(int)
    targets[3] -= 1
    on_targets_idx = targets[
        (targets[2] == targets['chrom']) & (targets[3] == targets['start'])
    ].index
    targets = targets.drop(on_targets_idx)

    # extract the number of mismatches in the off-target
    targets['mm'] = targets.apply(lambda x: x[13].split(':')[-1], axis=1)

    # count off-targets for a given guide
    mismatches_count = (
        targets.groupby(['id'])['mm'].value_counts().unstack().fillna(0).reset_index()
    )

    for m in ['0', '1', '2', '3']:
        if m not in mismatches_count:
            mismatches_count[m] = 0

    mismatches_count['id'] = mismatches_count['id'].astype(int)

    # print(mismatches_count)
    # add mismatch info to the guide dict
    for ix, x in mismatches_count.iterrows():
        i = int(x.id)
        counts = x[['0', '1', '2', '3']].to_dict()
        
        score_matrix = np.array([5, 3, 2, 1])
        count_values = np.array(list(counts.values()))

        cut_sites[i]['offtargets_sum_score'] = np.sum(score_matrix * count_values)
        cut_sites[i]['mm'] = counts
        cut_sites[i]['offtargets_str'] = '{:0.0f}|{:0.0f}|{:0.0f}|{:0.0f}'.format(
            *counts.values()
        )
        cut_sites[i]['offtargets_n'] = sum(counts.values())

        

    return cut_sites, targets
