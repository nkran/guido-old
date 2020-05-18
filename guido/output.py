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

def prepare_annotations(cut_sites):
    
    updated_cut_sites = []

    for cut_site in cut_sites:
        annotation_strings = []
        for ix, a in cut_site['annotation'].iterrows():
            if a['Feature'] == 'exon':
                label = '{}-E{}'.format(a['transcript_id'], a['Exon'])
            elif a['Feature'] in ['transcript', 'exon', 'CDS', 'start_codon', 'stop_codon', 'five_prime_utr', 'three_prime_utr']:
                label = '{}-{}'.format(a['Feature'], a['transcript_id'])
            else:
                label = '{}-{}'.format(a['Feature'], a['ID'])
            annotation_strings.append(label)

        cut_site['annotation_string'] = ' '.join(annotation_strings)
        cut_site['annotation_strings'] = annotation_strings

        updated_cut_sites.append(cut_site)

    return updated_cut_sites


def render_output(cut_sites, output_folder, targets_df=None):

    cut_sites = prepare_annotations(cut_sites)
    cs_df = pd.DataFrame(cut_sites).sort_values('absolute_cut_pos')

    cs_df[['chrom', 'start', 'end']] = pd.DataFrame(cs_df['guide_loc'].values.tolist())
    cs_df['variants'] = pd.DataFrame(cs_df['variants'])

    output = cs_df[
        [
            'chrom',
            'start',
            'end',
            'strand',
            'absolute_cut_pos',
            'guide',
            'cons_score',
            'variants_n',
            'complete_score',
            'sum_score',
            'offtargets_str',
            'offtargets_n',
            'annotation_string'
        ]
    ]
    output = output.sort_values(by=['chrom', 'start'])
    output.to_csv(os.path.join(output_folder, 'guides_list.csv'), index=False)

    if targets_df is not None:
        targets_df['diff'] = targets_df.apply(
            lambda x: parse_MD_tag(x[9], x[12]), axis=1
        )
        targets_df['id'] = targets_df['id'].astype(int)
        targets = [
            d
            for ix, d in targets_df[
                ['id', 'diff', 2, 3, 7, 9, 12, 'chrom', 'start', 'seq', 'mm']
            ]
            .to_dict(orient='index')
            .items()
        ]

        targets_sorted = sorted(targets, key=lambda x: x['id'])
        targets_grp = {
            k: list(g) for k, g in itertools.groupby(targets_sorted, lambda x: x['id'])
        }
    else:
        targets_grp = False

    template = env.get_template('guide_details.j2')
    output_details = template.render(
        cut_sites=cut_sites, targets_grp=targets_grp, rev_comp=rev_comp
    )

    with open(os.path.join(output_folder, 'guides_list_detailed.txt'), 'w') as f:
        f.write(output_details)
