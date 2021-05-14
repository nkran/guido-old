import pickle
import pandas as pd
import numpy as np
import allel
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import rankdata

import matplotlib as mpl
import seaborn as sns
import matplotlib.pyplot as plt

import h5py
import zarr
from scipy import stats
import malariagen_data
from decipy import executors as exe # dependency

def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in seq[::-1]])

# ag3 = malariagen_data.Ag3("gs://vo_agam_release/")
data = zarr.open('/Users/nace/imperial/conservation/data/AgamP4_conservation.zarr', mode='r')
root = zarr.open('/Users/nace/imperial/conservation/data/AgamP4_conservation_separated.zarr', mode='r')
accessibility = h5py.File('/Users/nace/imperial/conservation/data/accessibility.h5', mode='r')
# samples = ag3.sample_metadata(sample_sets="v3")

# genes_of_interest = ['AGAP004050', 'AGAP005958', 'AGAP007280', 'AGAP011377', 'AGAP013051']
# genes_of_interest = ['AGAP007165', 'AGAP006187', 'AGAP004203', 'AGAP000427', 'AGAP029113']
genes_of_interest = ['AGAP011515', 'AGAP008433', 'AGAP006385', 'AGAP005310']

def guide_accessibility(g, acc):
    chromosome, start, end = g['guide_loc']

    is_accessible = all(acc[chromosome]['is_accessible'][start-1:end])
    return is_accessible


def pam_snp_alleles(g):
    strand = g['strand']
    if strand == '-':
        pam_snp = [v['ac'] if v['pos_rel'] < 2 else 0 for v in g['variants']]
    else:
        pam_snp = [v['ac'] if v['pos_rel']>len(v['guide_seq'])-2 else 0 for v in g['variants']]

    return sum(pam_snp)


def seed_region_snp_alleles(g):
    strand = g['strand']
    if strand == '-':
        seed_snp = [v['ac'] if v['pos_rel'] > 2 and v['pos_rel'] < 12 else 0 for v in g['variants']]
    else:
        seed_snp = [v['ac'] if v['pos_rel']>len(v['guide_seq']) - 12 and v['pos_rel'] < len(v['guide_seq']) - 2 else 0 for v in g['variants']]

    return sum(seed_snp)


def small_seed_region_snp_alleles(g):
    strand = g['strand']
    if strand == '-':
        seed_snp = [v['ac'] if v['pos_rel'] > 2 and v['pos_rel'] < 6 else 0 for v in g['variants']]
    else:
        seed_snp = [v['ac'] if v['pos_rel']>len(v['guide_seq']) - 6 and v['pos_rel'] < len(v['guide_seq']) - 2 else 0 for v in g['variants']]

    return sum(seed_snp)


def parse_oof_deletions(mmej_patterns):
    return '|'.join([p for p in mmej_patterns['frame_shift']])


def get_ranking(G):
    ranking_columns = ['sum_score', 'complete_score', 'cons_score', 'offtargets_sum_score', 'variants_sum_allele_count', 
                       'variants_n', 'pam_snp_alleles', 'seed_region_snp_alleles', 'small_seed_region_snp_alleles']

    ranking_columns_beneficial = [True, True, True, False, False, False, False, False, False]
    weights = [1, 1, 1, 1, 1, 1, 1, 1, 1]
    # weights = [0.1, 0.2, 0.4, 0.5, 0.5, 0.5, 0.8, 0.4, 0.6]
    xij = G.loc[G.is_accessible,ranking_columns]

    kwargs = {
    'data': xij,
    'beneficial': ranking_columns_beneficial,
    'weights': weights,
    'rank_reverse': True,
    'rank_method': "ordinal"
    }

    # Build MCDM Executor
    wsm = exe.WSM(**kwargs) # Weighted Sum Method
    topsis = exe.Topsis(**kwargs) # Topsis 
    vikor = exe.Vikor(**kwargs) # Vikor

    analizer = exe.RankSimilarityAnalyzer()

    # Add MCDMs to anlizer
    analizer.add_executor(wsm)
    analizer.add_executor(topsis)
    analizer.add_executor(vikor)

    # run analizer
    G_ranks_weighted = analizer.get_ranks_dataframe()
    df = G.merge(G_ranks_weighted, how='left', left_index=True, right_index=True)

    return df


for gene_name in genes_of_interest:
    print(gene_name)
    G = pd.read_hdf(f'results/{gene_name}/guides_data.h5')
    # G = pd.read_hdf(f'results/{gene_name}/guides_data_{gene_name}.h5')
    G[['chrom', 'start', 'end']] = pd.DataFrame(G['guide_loc'].values.tolist())
    G['is_accessible'] = G.apply(lambda x: guide_accessibility(x, accessibility), axis=1)
    G['pam_snp_alleles'] = G.apply(lambda x: pam_snp_alleles(x), axis=1)
    G['has_pam_snp'] = G.apply(lambda x: x['pam_snp_alleles'] > 0, axis=1)
    G['seed_region_snp_alleles'] = G.apply(lambda x: seed_region_snp_alleles(x), axis=1)
    G['small_seed_region_snp_alleles'] = G.apply(lambda x: small_seed_region_snp_alleles(x), axis=1)
    G['mmej_oof'] = G.apply(lambda x: parse_oof_deletions(x.mmej_patterns), axis=1)
    G['variants_sum_allele_count'] = G.apply(lambda x: sum([v['ac'] for v in x['variants']]), axis=1)
    
    """ ------------ Ranking -------------- """
    G = get_ranking(G)
    
    G.to_hdf(f'results/{gene_name}/guides_data_{gene_name}.v5.h5', key='df', mode='w')
    G[
            [
                'chrom',
                'start',
                'end',
                'strand',
                'absolute_cut_pos',
                'guide',
                'cons_score',
                'variants_n',
                'variants_sum_allele_count',
                'pam_snp_alleles', 
                'seed_region_snp_alleles',
                'small_seed_region_snp_alleles',
                'sum_score',
                'complete_score',
                'mmej_oof',
                'offtargets_str',
                'offtargets_n',
                'offtargets_sum_score',
                'is_accessible',
                'has_pam_snp',
                'WSM',
                'Topsis',
                'Vikor',
                'annotation_string'
            ]
        ].to_csv(f'results/{gene_name}/guides_list_{gene_name}.v5.txt', sep='\t', index=True)


    # variation file
    def output_variation(variation):
        seq = variation[0]['guide_seq']
        out = '{}\t{}\t{}\n'.format(' '.join(variation[0]['guide_seq']), 'Freq', 'Allele count')
        for v in variation:
            s = '- '*(v['pos_rel']-1) + v['alt'] + ' ' + '- '*(len(seq) - v['pos_rel'])
            out += '{}\t{}\t\t{}\n'.format(s, round(v['af'],5), v['ac'])
        return out


    def output_subpop_variation(variation):
        seq = variation[0]['guide_seq']
        out = '{}\t{}\t{}\t{}\n'.format(' '.join(variation[0]['guide_seq']), 'Freq', 'Allele count', 'Country and species')

        for v in sorted(variation, key=lambda x: (x['pos_rel'], x['country'])):
            s = '- '*(v['pos_rel']-1) + v['alt'] + ' ' + '- '*(len(seq) - v['pos_rel'])
            species_string = ' '.join([f'{spc_tpl[0]}: {spc_tpl[1]*100:.1f}' for spc_tpl in v['species'] if spc_tpl[1] > 0.0])

            out += '{}\t{}\t{}\t{} ({})\n'.format(s, round(v['af'],5), v['ac'], v['country'], species_string)
        return out

    with open(f'results/{gene_name}/guides_ag3_variation_{gene_name}.v4.txt', 'w') as f:
        for ix, g in G.iterrows():
            f.write('gRNA-{}: {}\nLocation: {}-{}:{} ({})\n\n'.format(ix, g['guide'], g['guide_loc'][0], g['guide_loc'][1], g['guide_loc'][2], g['strand']))
            if len(g['variants']) > 0:
                f.write('Variation\n' + output_variation(g['variants']) + '---------\n')
                f.write('Regional variation\n' + output_subpop_variation(g['variants_subpop']) + '---------\n\n\n')
            else:
                f.write('No variation \n\n')
