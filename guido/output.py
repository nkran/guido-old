import os
import csv
import itertools
import pandas as pd

from math import trunc
from jinja2 import Environment, FileSystemLoader

from guido.helpers import parse_MD_tag, rev_comp
import guido.log as log

logger = log.createCustomLogger('output')
ROOT_PATH = os.path.abspath(os.path.dirname(__file__))

file_loader = FileSystemLoader(os.path.join(ROOT_PATH, 'templates'))
env = Environment(loader=file_loader)
env.filters['rev_comp'] = rev_comp


def render_output(cut_sites, targets_df, output_folder):

    cs_df = pd.DataFrame(cut_sites)

    cs_df[['chrom', 'start', 'end']] = pd.DataFrame(cs_df['guide_loc'].values.tolist())
    cs_df['variants'] = pd.DataFrame(cs_df['variants'])

    output = cs_df[['chrom', 'start', 'end', 'strand', 'absolute_cut_pos', 'guide', 'cons_score', 'variants_n', 'complete_score', 'sum_score', 'offtargets_str', 'offtargets_n']]
    output = output.sort_values(by=['chrom', 'start'])
    output.to_csv(os.path.join(output_folder, 'guides_list.csv'), index=False)

    targets_df['diff'] = targets_df.apply(lambda x: parse_MD_tag(x[9], x[12]), axis=1)
    targets_df['id'] = targets_df['id'].astype(int)
    targets = [d for ix, d in targets_df[['id', 'diff', 2, 3, 7, 9, 12, 'chrom', 'start', 'seq', 'mm']].to_dict(orient='index').items()]
    
    targets_sorted = sorted(targets, key=lambda x: x['id'])
    targets_grp = {k: list(g) for k, g in itertools.groupby(targets_sorted, lambda x: x['id'])}
    template = env.get_template('guide_details.j2')
    output_details = template.render(cut_sites=cut_sites, targets_grp=targets_grp, rev_comp=rev_comp)

    with open(os.path.join(output_folder, 'guides_list_detailed.txt'), 'w') as f:
        f.write(output_details)

