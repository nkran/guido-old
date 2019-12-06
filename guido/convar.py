import os
import subprocess
import h5py
import zarr
import allel
import numpy as np
import pandas as pd

import guido.log as log

logger = log.createCustomLogger('Conservation and variation')
ROOT_PATH = os.path.abspath(os.path.dirname(__file__))


def fetch_convar_score(cut_site, conservation, var_callset):
    '''
    Returns mean conservation score and evaluates SNPs within the guide site
    '''

    # Calculate conservation score
    chrom, start, end = cut_site['guide_loc']

    if conservation:
        cut_site['cons_score'] = np.mean(conservation['joined/{}/score/'.format(chrom)][3,start:end])

    if var_callset:
        variants = var_callset[chrom]['variants']
        calldata = var_callset[chrom]['calldata']
        pos = allel.SortedIndex(variants['POS'])

        cut_site['variants'] = {}
        cut_site['variants_n'] = 0

        try:
            loc = pos.locate_range(start, end)
            g = allel.GenotypeArray(calldata['GT'][loc])

            cut_site['variants'] = {
                'pos': variants['POS'][loc],
                'alt': variants['ALT'][loc],
                'ref': variants['REF'][loc],
                'allel_count': np.array(g.count_alleles())
                }
            cut_site['variants_zipped'] = list(zip(*cut_site['variants'].values()))
            cut_site['variants_n'] = len(variants['POS'][loc])

        except:
            pass

    return cut_site


def apply_conservation_variation_score(cut_sites, conservation_store, variation_store, pool):
    '''
    STORE STRUCTURE:
    scores
        - 2L
            - AaraD1
                - identity (8, 2L len)
                    - '27_30'
                    - '29_30'
                    - '30_30'
                    - '35_50'
                    - '45_50'
                    - '48_50'
                    - '49_50'
                    - '50_50'
                - score (2, 2L len)
                    - short window average (mean of 27, 29, 30 / 30)
                    - long window average (mean of 45, 48, 49, 50 / 50)
                    - average
    joined
        - 2L
            - averages (21, 2L len)
                - AaegL5
                - AalbS2
                - AaraD1
                ...
            - weighted_averages (21, 2L len)
                - ... each genome
            - score (5, 2L len)
                - average of short + long across all species - weighted by *above
                - SNPs
                - varscore - conservation score combined with SNP density
                - varscore scaled to 0-1 scale
                - varscore non-coding (TODO)

            - rolling (8, 2L len) - calculate on the go?
                - SNP density
                - short window average - rolling (20 bp)
                - long window average - rolling (5 kb)
                - average of short + long across all species - rolling (20 bp)
                - average of short + long across all species - rolling (5 kb)
                - varscore - rolling (20 bp)
                - varscore - rolling (5 kb)
                - varscore non-coding - rolling (20b)
                - varscore non-coding - rolling (5 kb)

        --- ### calculate rolling means ### ---
            mean_short_sliding = np.array(pd.DataFrame(np.mean(identity_short, axis=0)).rolling(20, win_type ='hamming', center=True).mean())
            mean_long_sliding = np.array(pd.DataFrame(np.mean(identity_long, axis=0)).rolling(20, win_type ='hamming', center=True).mean())

            mean_short_sliding_big = np.array(pd.DataFrame(np.mean(identity_short, axis=0)).rolling(5000, win_type ='hamming', center=True).mean())
            mean_long_sliding_big = np.array(pd.DataFrame(np.mean(identity_long, axis=0)).rolling(5000, win_type ='hamming', center=True).mean())
    '''

    '''
        2L
            calldata
                GT (1497107, 20, 2) int8
            samples (20,) object
            variants
                ALT (1497107, 3) object
                CHROM (1497107,) object
                FILTER_PASS (1497107,) bool
                ID (1497107,) object
                POS (1497107,) int32
                QUAL (1497107,) float32
                REF (1497107,) object
    '''

    # handle variation file type
    if variation_store is not False and os.path.exists(variation_store):
        var_path, var_ext = os.path.splitext(variation_store)

        if var_ext == '.vcf':
            logger.info('Reading VCF file. To speed up the process transform your file into HDF5 or Zarr format.')
            var_callset = allel.read_vcf(variation_store)
        elif var_ext == '.h5' or var_ext == '.hdf5':
            var_callset = h5py.File(variation_store, mode='r')
        elif var_ext == '.zarr':
            var_callset = zarr.open(variation_store, 'r')
        else:
            logger.error('Variation file could not be opened. Please check the supported file formats.')
            quit()
    else:
        var_callset = None
        logger.info('Variation file was not provided.')

    # open conservation store
    if conservation_store is not False and os.path.exists(conservation_store):
        conservation = zarr.open(conservation_store, 'r')
    else:
        conservation = None
        logger.info('Conservation resource was not provided.')

    if any([var_callset, conservation]):
        iterable_cut_sites = ((cut_site, conservation, var_callset) for cut_site in cut_sites)
        cut_sites = pool.starmap(fetch_convar_score, iterable_cut_sites)

    return cut_sites