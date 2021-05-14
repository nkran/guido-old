import numpy as np
import allel

import malariagen_data
import dask 

ag3 = malariagen_data.Ag3("gs://vo_agam_release/")
samples = ag3.sample_metadata(sample_sets="v3")

guides_var = []
country_samples = {l: list(cdf.index) for l, cdf in samples.groupby('country')}

df = ag3.geneset()

genes = ['AGAP011515', 'AGAP008433']
genes = ['AGAP011515', 'AGAP008433', 'AGAP006385', 'AGAP005310']

for gid in genes:
    print(gid)
    gene = df.loc[df.ID == gid,:]
    chromosome = gene['seqid'].to_string(index=False).strip()
    pos, ref, alt = ag3.snp_sites(contig=chromosome, site_mask="gamb_colu_arab")
    gt = ag3.snp_genotypes(contig=chromosome, sample_sets="v3", site_mask="gamb_colu_arab")
    gt = allel.GenotypeDaskArray(gt)
    pos = allel.SortedIndex(pos)

    loc_gene = pos.locate_range(int(gene.start), int(gene.end))
    loc_pos = pos[loc_gene]

    with open(f'pos_{gid}.npy', 'wb') as f:
        np.save(f, pos[loc_gene].reshape(1, -1))

    with open(f'ref_{gid}.npy', 'wb') as f:
        np.save(f, ref[loc_gene].compute().reshape(1, -1))
    
    with open(f'alt_{gid}.npy', 'wb') as f:
        np.save(f, alt[loc_gene].compute())
    
    with open(f'loc_{gid}.npy', 'wb') as f:
        np.save(f, loc_gene)

    gt[loc_gene].rechunk().to_zarr(f'_gt_{gid}.zarr', overwrite=True)
