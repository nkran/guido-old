from IPython import get_ipython

import pickle
import pandas as pd
import numpy as np
import allel

import matplotlib as mpl
import seaborn as sns
import matplotlib.pyplot as plt

import h5py
import zarr
from scipy import stats
import malariagen_data
import humanize

def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in seq[::-1]])

def output_variation(variation):
    seq = variation[0]['guide_seq']

    if variation[0]['guide_strand'] == '-':
        ref = ' '.join(rev_comp(variation[0]['guide_seq']))
    else:
        ref = ' '.join(variation[0]['guide_seq'])
    out = '{}\t{}\t{}\n'.format(ref, 'Freq', 'Allele count')
    for v in variation:
        s = '- '*(v['pos_rel']-1) + v['alt'] + ' ' + '- '*(len(seq) - v['pos_rel'])
        out += '{}\t{}\t\t{}\n'.format(s, round(v['af'],5), v['ac'])
    return out


def output_subpop_variation(variation):
    seq = variation[0]['guide_seq']

    if variation[0]['guide_strand'] == '-':
        ref = ' '.join(rev_comp(variation[0]['guide_seq']))
    else:
        ref = ' '.join(variation[0]['guide_seq'])
    out = '{}\t{}\t{}\t{}\n'.format(ref, 'Freq', 'Allele count', 'Country and species')

    for v in sorted(variation, key=lambda x: (x['pos_rel'], x['country'])):
        s = '- '*(v['pos_rel']-1) + v['alt'] + ' ' + '- '*(len(seq) - v['pos_rel'])
        species_string = ' '.join([f'{spc_tpl[0]}: {spc_tpl[1]*100:.1f}' for spc_tpl in v['species'] if spc_tpl[1] > 0.0])

        out += '{}\t{}\t{}\t{} ({})\n'.format(s, round(v['af'],5), v['ac'], v['country'], species_string)
    return out


def map_alt_alleles(guide_pos, guide_pos_rel, guide_alt, guide_ac, guide_dict, country=None, species=None):
        
        variation = []
        guide_af = guide_ac.to_frequencies()

        for iy, alt in enumerate(guide_alt):
            for ix, af in enumerate(guide_af[iy,1:]):
                if float(af) > 0:
                    v = {}
                    v['pos_rel'] = guide_pos_rel[iy]
                    v['pos'] = guide_pos[iy]
                    v['ref'] = guide_ref[iy].decode("utf-8")
                    v['alt'] = alt[ix].decode("utf-8")
                    v['ac'] = guide_ac[iy,ix+1]
                    v['af'] = float(af)
                    v['aft'] = guide_ac[iy,ix+1] / (3081*2)
                    v['country'] = country

                    if type(species) == dict:
                        total = guide_ac[iy,ix+1]
                        v['species'] = [(s, counts[iy,ix+1] / total) for s, counts in species.items()]
                    else:
                        v['species'] = species
                        
                    v['guide_ix'] = df.name

                    if guide_dict['strand'] == '-':
                        v['guide_seq'] = rev_comp(guide_dict['guide'])
                    else:
                        v['guide_seq'] = guide_dict['guide']

                    v['guide_strand'] = guide_dict['strand']
                    variation.append(v)

        return variation


def plot_transcripts(geneset, chrom, start, stop, height=.2, label_transcripts=True, label_exons=False, label_exon_size=False,
                     label_codons=False, highlight_exons=None, label_cdss=False, highlight_color='red', ax=None,
                     title=None):
    """Plot all transcripts for all genes overlapping a given region.
        Authored by: Alistar Miles (https://alimanfoo.github.io/) """

    if ax is None:
        fig, ax = plt.subplots(figsize=(14, 6))
        sns.despine(ax=ax, left=True, offset=5)
        
    if title:
        ax.set_title(title, va='bottom')

    # find genes overlapping the given region 
    genes = geneset.query("((type == 'gene') or (type == 'ncRNA_gene')) and (seqid == %r) and (end >= %s) and (start <= %s)" % (chrom, start, stop)).sort_values('start')

    # iterate over genes
    for _, gene in genes.iterrows():

        # find child transcripts
        transcripts = geneset.query("(type == 'mRNA') and (Parent == %r)" % gene.ID).sort_values('ID')

        if gene.type == 'ncRNA_gene':
            nc_transcripts = geneset.query("(type == 'ncRNA_gene') and (seqid == %r) and (end >= %s) and (start <= %s)" % (chrom, start, stop)).sort_values('start')
            transcripts = pd.concat([nc_transcripts, transcripts]).sort_values('ID')
        
        # iterate over transcripts
        for i, (_, transcript) in enumerate(transcripts.iterrows()):
            # coordinates for plotting the transcript
            if transcript.strand == '+':
                y = i
            else:
                y = -i - 1

            # annotate with transcript ID
            text_y = y + height + (height / 10)
            if label_transcripts == 'right':
                text_x = min(stop, transcript.end)
                ha = 'right'
            else:
                text_x = max(start, transcript.start)
                ha = 'left'
            if label_transcripts:
                if transcript.strand == '+':
                    text = '%s >' % transcript.ID
                else:
                    text = '< %s' % transcript.ID
                ax.text(text_x, text_y, text, ha=ha, va='bottom')
                        
            # find child exons
            exons = geneset.query("type == 'exon' and Parent == %r" % transcript.ID).sort_values('start')

            # iterate over exons to plot introns
            last_exon = None
            for i, (_, exon) in enumerate(exons.iterrows()):
                x = exon.start
                width = exon.end - x
                # plot intron
                if last_exon is not None:
                    ax.plot([last_exon.end, (last_exon.end + exon.start) / 2, exon.start], [y + height / 2, y + height / 1.5, y + height / 2], 'k-')
                last_exon = exon
                
                # exon number
                n = i + 1 if exon.strand == '+' else len(exons) - i
                
                # label exons
                if label_exons and exon.end > start and exon.start < stop:
                    text_x = (exon.start + exon.end) / 2
                    ha = 'center'
                    if text_x < start:
                        text_x = start
                        ha = 'left'
                    elif text_x > stop:
                        text_x = stop
                        ha = 'right'
                    s = str(n)
                    if label_exon_size:
                        s += ' (%s)' % (exon.end - exon.start + 1)
                    if label_exons == 'center':
                        ax.text(text_x, y + height / 2, s, ha=ha, va='center', color='w', zorder=20, fontweight='bold')
                    else:
                        ax.text(text_x, text_y, s, ha=ha, va='bottom', color='k', zorder=20)
                
                # highlight exons
                if highlight_exons and (transcript.ID, n) in highlight_exons:
                    patch = plt.Rectangle((x, y), width, height, color=highlight_color, alpha=0.4, zorder=10)
                    ax.add_patch(patch)

            # find child CDSs
            cdss = geneset.query("(type == 'CDS' or type == 'rRNA') and Parent == %r" % transcript.ID)
            if transcript.strand == '+':
                cdss = cdss.sort_values('start', ascending=True)
            else:
                cdss = cdss.sort_values('end', ascending=False)
            
            # keep track of CDS position
            cds_pos = 0
            
            # plot CDSs
            for _, cds in cdss.iterrows():
                x = cds.start
                width = cds.end - x
                
                # plot CDS
                patch = plt.Rectangle((x, y), width, height, color='k')
                ax.add_patch(patch)
                
                if label_codons:
                    # report 1-based numbers
                    s = '%s (%s)' % ((cds_pos // 3) + 1, cds_pos + 1)
                    if transcript.strand == '+':
                        text_x = x
                        ha = 'left'
                    else:
                        text_x = x + width
                        ha = 'right'
                    if text_x > start and text_x < stop:
                        ax.text(text_x, text_y, s, ha=ha, va='bottom')
                                
                # label CDSs
                if label_cdss and cds.end > start and cds.start < stop:
                    text_x = (cds.start + cds.end) / 2
                    ha = 'center'
                    if text_x < start:
                        text_x = start
                        ha = 'left'
                    elif text_x > stop:
                        text_x = stop
                        ha = 'right'
                    s = cds.ID
                    if label_cdss == 'center':
                        ax.text(text_x, y + height / 2, s, ha=ha, va='center', color='w', zorder=20, fontweight='bold')
                    else:
                        ax.text(text_x, text_y, s, ha=ha, va='bottom', color='k', zorder=20)
                
                # accumulate CDS positions
                cds_pos += width + 1  # N.B., GFF coords are 1-based end-inclusive

            # find child UTRs
            utrs = geneset.query("(type == 'three_prime_UTR' or type == 'five_prime_UTR') and Parent == %r" % transcript.ID).sort_values('start')
            for _, utr in utrs.iterrows():
                x = utr.start
                width = utr.end - x
                utr_height = height * .8
                utr_y = y + (height - utr_height) / 2
                patch = plt.Rectangle((x, utr_y), width, utr_height, facecolor='#cccccc', edgecolor='k')
                ax.add_patch(patch)

    ax.set_yticks([])
    ax.set_xlim(start, stop)
    ax.set_xticklabels([humanize.intcomma(int(x)) for x in ax.get_xticks()])
    ax.axhline(0 - (height / 2), color='k', linestyle='--')
    ax.set_xlabel('Chromosome %s position (bp)' % chrom)
    ax.autoscale(axis='y', tight=True)
    
    return ax


def geneset_to_pandas(geneset):
    """Life is a bit easier when a geneset is a pandas DataFrame."""
    items = []
    for n in geneset.dtype.names:
        v = geneset[n]
        # convert bytes columns to unicode (which pandas then converts to object)
        if v.dtype.kind == 'S':
            v = v.astype('U')
        items.append((n, v))
    return pd.DataFrame.from_dict(dict(items))


def plot_guides(guides, start, ax):
    # guides = [(start, end, colour)]
    for guide in guides:
        if len(guide) == 3:
            colour = guide[2]
        else:
            colour = 'orange'

        ax.axvspan(guide[0] - start, guide[1] - start, ymax=95, color=colour, alpha=0.5)


def plot_gene_conservation(gene, output=None):
    df_gene = geneset_agam.loc[geneset_agam.ID == gene,:].to_dict(orient='records')
    x_chromosome = df_gene[0]['seqid'] 
    x_start = df_gene[0]['start'] - 100
    x_end = df_gene[0]['end'] + 100

    x_cs = data['{}/Cs'.format(x_chromosome)][0, x_start:x_end]
    
    # x_stack = data['{}/stack'.format(x_chromosome)][:,x_start:x_end]
    # x_stack_labels = data['{}/stack'.format(x_chromosome)].attrs['species']

    x_accessibility = accessibility['/{}/is_accessible'.format(x_chromosome)][x_start:x_end]
    x_acc = np.ones(len(x_cs))
    x_acc[x_accessibility] = 0

    fig, ax = plt.subplots(2, 1, figsize=(20, 7))
    plot_1 = ax[0]
    plot_2 = ax[1]
    # plot_3 = ax[2]
    # plot transcripts
    plot_transcripts(geneset_agam, x_chromosome, x_start, x_end, label_codons=False, label_exons=True, ax=plot_1)
    sns.despine(ax=plot_1, left=True, offset=5)  

    plot_2.plot(x_cs, color='gray', linewidth=2, alpha=1)
    plot_2.vlines(np.where(x_acc), ymin=0, ymax=1, color='gray', alpha=.05)

    plot_2.set_xlim(0, len(x_cs))
    plot_2.set_ylim(0,1)

    # plt.show()

    if output:
        fig.savefig(f'{output}/{gene}.conservation.pdf', dpi=300) 


plt.rcParams['figure.figsize'] = (25.0, 10.0)

identities =    ['27_30', '29_30', '30_30', '35_50', '45_50', '48_50', '49_50', '50_50']
genomes =       ['AaegL5', 'AalbS2', 'AaraD1', 'AatrE3', 'AchrA1', 'AcolM1', 'AculA1', 'AdarC3', 
                 'AdirW1', 'AepiE1', 'AfarF2', 'AfunF1', 'AmacM1', 'AmelC2', 'AmerM2', 'AminM1', 
                 'AquaS1', 'AsinC2', 'AsteI2', 'CpipJ2', 'DmelP6']
chromosomes =   ['2L', '2R', '3L', '3R', 'X']

# Load geneset and save it to pandas dataframe
geneset_agam = allel.FeatureTable.from_gff3('/Users/nace/imperial/conservation/data/ref/AgamP4.12.gff3',
                                            attributes=['ID', 'Parent', 'Name'])

geneset_agam = geneset_to_pandas(geneset_agam)

ag3 = malariagen_data.Ag3("gs://vo_agam_release/")
data = zarr.open('/Users/nace/imperial/conservation/data/AgamP4_conservation.zarr', mode='r')
root = zarr.open('/Users/nace/imperial/conservation/data/AgamP4_conservation_separated.zarr', mode='r')
accessibility = h5py.File('/Users/nace/imperial/conservation/data/accessibility.h5', mode='r')
samples = ag3.sample_metadata(sample_sets="v3")

# female fertility + sex diff
# genes_of_interest = ['AGAP004050', 'AGAP005958', 'AGAP007280', 'AGAP011377', 'AGAP013051', 'AGAP029421']
# genes_of_interest = ['AGAP004050', 'AGAP007280', 'AGAP011377', 'AGAP013051', 'AGAP029421']
# genes_of_interest = ['AGAP005958']

# nicole + will
# genes_of_interest = ['AGAP007165', 'AGAP006187', 'AGAP004203', 'AGAP000427', 'AGAP029113']
# genes_of_interest = [ 'AGAP000427', 'AGAP029113']

genes_of_interest = ['AGAP011515', 'AGAP008433', 'AGAP006385', 'AGAP005310']

for gene_name in genes_of_interest:
    print(gene_name, '------------------\n')
    with open(f'results/{gene_name}/cut_sites_all.pickle', 'rb') as p:
        guides = pickle.load(p)

    plot_gene_conservation(gene_name, output=f'/Users/nace/imperial/guides/results/{gene_name}/')

    G = pd.DataFrame(guides)
    G.loc[pd.isna(G['offtargets_sum_score']),'offtargets_sum_score'] = 0
    df = G.copy()
    # df.sort_values('cons_score', ascending=False, inplace=True)
    # df = df.loc[(df['cons_score'] > 0.3)&(df['offtargets_n'] < 3),:]

    gt = zarr.open(f'results/var/_gt_{gene_name}.zarr', 'r')
    alt = np.load(f'results/var/nps/alt_{gene_name}.npy')
    ref = np.load(f'results/var/nps/ref_{gene_name}.npy')
    pos = np.load(f'results/var/nps/pos_{gene_name}.npy')
    loc = np.load(f'results/var/nps/loc_{gene_name}.npy', allow_pickle=True)

    pos = allel.SortedIndex(pos[0])
    gt = allel.GenotypeChunkedArray(gt)

    country_samples = {l: list(cdf.index) for l, cdf in samples.groupby('country')}
    variation = []

    for ix, df in G.iterrows():
        
        guide = df.to_dict()
        start, end = guide['guide_loc'][1], guide['guide_loc'][2]

        try:
            guide_loc = pos.locate_range(start + 1, end)
            guide_pos = pos[guide_loc]
            guide_alt = alt[guide_loc]
            guide_ref = ref[0][guide_loc]
            guide_gt = gt[guide_loc]

            is_variant = guide_gt.count_alleles().is_variant()
            # select only loci with alt variants
            guide_pos = guide_pos[is_variant]

            if len(guide_pos) > 0:
                guide_pos_rel = guide_pos - start
                guide_pos_rel = allel.SortedIndex(guide_pos_rel)
                guide_alt = guide_alt[is_variant,:]
                guide_ref = guide_ref[is_variant]

                guide_ac = guide_gt[is_variant,:,:].count_alleles()

                v = map_alt_alleles(guide_pos, guide_pos_rel, guide_alt, guide_ac, df)
                variation.extend(v)

        except KeyError:
            print('No variation for ', ix, start, end)
        

    vdf = pd.DataFrame(variation)
    v_grouped = []

    for country, samples_df in samples.groupby('country'):
        
        print(country, len(samples_df))
        group_ix = samples_df.index.values
        species_group_ix = {}

        for s in set(samples_df['species']):
            species_group_ix[s] = samples_df.query(f'species == "{s}"').index.values
        
        for ix, df in G.iterrows():
            
            guide = df.to_dict()
            start, end = guide['guide_loc'][1], guide['guide_loc'][2]

            try:
                guide_loc = pos.locate_range(start + 1, end)
                guide_pos = pos[guide_loc]
                guide_alt = alt[guide_loc]
                guide_ref = ref[0][guide_loc]
                guide_gt = allel.GenotypeDaskArray(gt)[guide_loc, group_ix].compute()
                guide_gt = allel.GenotypeChunkedArray(guide_gt)

                guide_gt_subpops = {subpop: allel.GenotypeChunkedArray(allel.GenotypeDaskArray(gt)[guide_loc, ix].compute()) for subpop, ix in species_group_ix.items()}

                is_variant = guide_gt.count_alleles(max_allele=3).is_variant()
                # select only loci with alt variants
                guide_pos = guide_pos[is_variant]
                guide_pos_rel = guide_pos - start
                guide_pos_rel = allel.SortedIndex(guide_pos_rel)

                guide_alt = guide_alt[is_variant,:]
                guide_ref = guide_ref[is_variant]

                    # species_allel_counts[s] = allel.GenotypeDaskArray(guide_gt)[:,species_ix]

                if len(guide_pos) > 0:
                    guide_ac = guide_gt[is_variant,:,:].count_alleles(max_allele=3)
                    guide_ac_subpop = {subpop: gt[is_variant,:,:].count_alleles(max_allele=3) for subpop, gt in guide_gt_subpops.items()}

                    v = map_alt_alleles(guide_pos, guide_pos_rel, guide_alt, guide_ac, df, country, guide_ac_subpop)
                    v_grouped.extend(v)

            except KeyError:
                pass
                # print('No variation for ', ix, start, end)


    vg_df = pd.DataFrame(v_grouped)

    G['variants'] = G.apply(lambda x: vdf.loc[vdf.guide_ix == x.name,:].to_dict(orient='records'), axis=1)
    G['variants_n'] = G.apply(lambda x: len(x['variants']), axis=1)
    G['variants_subpop'] = G.apply(lambda x: vg_df.loc[vg_df.guide_ix == x.name,:].to_dict(orient='records'), axis=1)

    with open(f'results/{gene_name}/guides_ag3_variation.txt', 'w') as f:
        for ix, g in G.iterrows():
            f.write('Guide: {}\nLocation: {}-{}:{} ({})\n\n'.format(g['guide'], g['guide_loc'][0], g['guide_loc'][1], g['guide_loc'][2], g['strand']))
            if len(g['variants']) > 0:
                f.write('Variation\n' + output_variation(g['variants']) + '---------\n')
                f.write('Regional variation\n' + output_subpop_variation(g['variants_subpop']) + '---------\n\n\n')
            else:
                f.write('No variation \n\n')



    G.to_pickle(f'results/{gene_name}/guides_ag3_variation.pickle')
    G.to_hdf(f'results/{gene_name}/guides_data.h5', key='df', mode='w')
    vg_df.to_hdf(f'results/{gene_name}/variation_grouped.h5', key='df', mode='w')
    vdf.to_hdf(f'results/{gene_name}/variation_overall.h5', key='df', mode='w')