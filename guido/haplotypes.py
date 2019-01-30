from __future__ import print_function

import h5py
import zarr
import numpy as np
import pandas as pd
import allel

from sys import stdout

# callset = zarr.open('/data2/archive/ag1000g/haps/ag1000g.phase2.ar1.haplotypes.2L.zarr', mode='r')

# gt_zarr = callset['2L/calldata/genotype']

# pos = allel.SortedIndex(callset['2L/variants/POS'])
# loc_region = pos.locate_range(44994010, 44994033)

# gt_region = allel.GenotypeArray(gt_zarr[loc_region])
# ref_region = callset['2L/variants/REF'][loc_region]
# alt_region = callset['2L/variants/ALT'][loc_region]
# pos_region = callset['2L/variants/POS'][loc_region]

# h = gt_region.to_haplotypes()
# n_snps, n_haps = h.shape

# h_df = pd.DataFrame(dict(zip(pos_region, h)))

# print ref_region
# print alt_region



# hap_count = h_df.groupby(h_df.columns.tolist(), as_index=True).size()

# hap_count_d = dict(hap_count)

# for hap, count in hap_count_d.items():
#     for i, allele in enumerate(hap):
#         if allele == 0:
#             print ref_region[i],
#         elif allele == 1:
#             print alt_region[i],
    
#     print count, round(float(count) / float(n_haps), 4), '----'

# samples_meta = pd.read_csv('/data2/archive/ag1000g/haps/haplotypes.meta.txt', sep='\t')
# samples_pops = samples_meta[['ox_code', 'population']].drop_duplicates()

# samples_callset = pd.DataFrame(zip(range(1, 2328), callset['2L/samples']))
# samples_callset = samples_callset.reindex(np.repeat(samples.index.values, 2)).reset_index(drop=True)
# samples_callset = samples_callset.dropna()

# h_samples_df = pd.concat([h_df, samples_callset], axis=1, sort=False)
# h_samples_df = pd.merge(h_samples_df, samples_pops, how='inner', left_on=1, right_on='ox_code')

# h_samples_df['hap'] = h_samples_df.apply(lambda x: ''.join([str(a) for a in x[0:len(pos_region)]]), axis=1)

# haps_sizes = h_samples_df.groupby(by='hap').size()
# haps_samples = h_samples_df.groupby(by='hap')['population'].agg(list)

# from collections import Counter
# hap_populations = [(hap, dict(Counter(pops))) for hap, pops in haps_samples.items()]

# for hap in hap_populations:
#     print '-----'
#     print hap[0]    
#     print '-'
#     for pop, n in hap[1].items():
#         print pop, n, round(float(n) / (pop_sizes.loc[pop]['count'] * 2), 3)

