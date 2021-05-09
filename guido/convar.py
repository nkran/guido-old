import os
import h5py
import zarr
import allel
import numpy as np

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
        cut_site['cons_score'] = np.mean(
            # root[chromosome][args.array]
            conservation[f'{chrom}/Cs'][0,start:end]
        )

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
                'allel_count': np.array(g.count_alleles()),
            }
            cut_site['variants_zipped'] = list(zip(*cut_site['variants'].values()))
            cut_site['variants_n'] = len(variants['POS'][loc])

        except Exception:
            pass

    return cut_site


def apply_conservation_variation_score(
    cut_sites, conservation_store, variation_store, pool
):

    # handle variation file type
    if variation_store is not False and os.path.exists(variation_store):
        var_path, var_ext = os.path.splitext(variation_store)

        if var_ext == '.vcf':
            logger.info(
                'Reading VCF file. To speed up the process transform your file into HDF5 or Zarr format.'
            )
            var_callset = allel.read_vcf(variation_store)
        elif var_ext == '.h5' or var_ext == '.hdf5':
            var_callset = h5py.File(variation_store, mode='r')
        elif var_ext == '.zarr':
            var_callset = zarr.open(variation_store, 'r')
        else:
            logger.error(
                'Variation file could not be opened. Please check the supported file formats.'
            )
            quit()
    else:
        var_callset = None
        logger.info('Variation file was not provided.')

    # open conservation store
    if conservation_store is not False and os.path.exists(conservation_store):
        ar_path, con_ext = os.path.splitext(conservation_store)

        if con_ext == '.zarr':
            conservation = zarr.open(conservation_store, 'r')
        elif con_ext == '.h5' or con_ext == '.hdf5':
            conservation = h5py.File(conservation_store, mode='r')
    else:
        conservation = None
        logger.info('Conservation resource was not provided.')

    if any([var_callset, conservation]) and con_ext == '.zarr':
        iterable_cut_sites = (
            (cut_site, conservation, var_callset) for cut_site in cut_sites
        )
        cut_sites = pool.starmap(fetch_convar_score, iterable_cut_sites)
    else:
        cut_sites = [fetch_convar_score(cut_site, conservation, var_callset) for cut_site in cut_sites]

    return cut_sites
