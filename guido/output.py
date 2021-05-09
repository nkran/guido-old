import os
import csv
import itertools
import pandas as pd
import pickle

from math import trunc
from jinja2 import Environment, FileSystemLoader

from guido.helpers import parse_MD_tag, rev_comp, parse_oof_deletions
from guido.ranking import rank_guides
import guido.log as log

logger = log.createCustomLogger('output')
ROOT_PATH = os.path.abspath(os.path.dirname(__file__))

file_loader = FileSystemLoader(os.path.join(ROOT_PATH, 'templates'))
env = Environment(loader=file_loader)
env.filters['rev_comp'] = rev_comp

def prepare_annotations(cut_sites, ann_ext):
    updated_cut_sites = []
    feature_dict = [
        'transcript',
        'CDS',
        'start_codon',
        'stop_codon',
        'five_prime_utr',
        'five_prime_UTR',
        'three_prime_utr',
        'three_prime_UTR'
    ]

    for cut_site in cut_sites:
        annotation_strings = []
        for ix, a in cut_site['annotation'].iterrows():
            if ann_ext in ['.gff3']:
                if a['Feature'] == 'gene':
                    label = 'gene-{}'.format(a['ID'])
                elif a['Feature'] == 'exon':
                    label = a['Exon']
                elif a['Feature'] == 'mRNA':
                    label = 'transcript-{}'.format(a['ID'])
                elif a['Feature'] in feature_dict:
                    label = '{}-{}'.format(a['Feature'], a['Parent'])
            elif ann_ext in ['.gtf']:
                if a['Feature'] == 'gene':
                    label = 'gene-{}'.format(a['ID'])
                elif a['Feature'] == 'exon':
                    label = '{}-E{}'.format(a['transcript_id'], a['Exon'])
                elif a['Feature'] in feature_dict:
                    label = '{}-{}'.format(a['Feature'], a['transcript_id'])
            annotation_strings.append(label)

        # gff3 file has every annotation twice -> consider using set(annotation_strings) although it loses previous intuitive order
        cut_site['annotation_string'] = ' '.join(annotation_strings) 
        cut_site['annotation_strings'] = annotation_strings

        updated_cut_sites.append(cut_site)

    return updated_cut_sites


def render_output(cut_sites, output_folder, ann_ext, targets_df=None):

    if ann_ext is not None:
        cut_sites = prepare_annotations(cut_sites, ann_ext)
    cs_df = pd.DataFrame(cut_sites).sort_values('absolute_cut_pos')

    cs_df[['chrom', 'start', 'end']] = pd.DataFrame(cs_df['guide_loc'].values.tolist())
    cs_df['variants'] = pd.DataFrame(cs_df['variants'])
    try:
        cs_df['mmej_oof'] = cs_df.apply(lambda x: parse_oof_deletions(x.mmej_patterns), axis=1)
    except:
        cs_df['mmej_oof'] = None

    cs_df = rank_guides(cs_df)

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
            'sum_score',
            'complete_score',
            'mmej_oof',
            'offtargets_str',
            'offtargets_n',
            'offtargets_sum_score',
            'weighted_sum',
            'rank'
        ]
    ]
    if ann_ext is not None:
        output.insert(len(output.columns), 'annotation_string', cs_df['annotation_string'])
    output = output.sort_values(by=['chrom', 'start'])
    output.to_csv(os.path.join(output_folder, 'guides_list.csv'), sep='\t', index=False)

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

    # with open(os.path.join(output_folder, 'targets_all.pickle'), 'wb') as p:
    #     pickle.dump(targets_grp, p)

    with open(os.path.join(output_folder, 'cut_sites_all.pickle'), 'wb') as p:
        pickle.dump(cut_sites, p)

    with open(os.path.join(output_folder, 'guides_list_detailed.txt'), 'w') as f:
        f.write(output_details)
