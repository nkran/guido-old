import os
import subprocess
import zarr
import numpy as np
import pandas as pd
from tqdm import tqdm
import guido.log as log

ROOT_PATH = os.path.abspath(os.path.dirname(__file__))


def fetch_convar_score(cut_site, z, variation):
    '''
    Returns mean conservation score and evaluates SNPs within the guide site
    '''

    # TODO: add positional variation

    chrom, start, end = cut_site['guide_loc']
    cut_site['cons_score'] = np.mean(z['joined/{}/score/'.format(chrom)][3,start:end])

    return cut_site


def apply_conservation_variation_score(cut_sites, store_location, variation, pool):
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
                - AatrE3 
                - AchrA1 
                - AcolM1 
                - AculA1 
                - AdarC3 
                - AdirW1
                - AepiE1
                - AfarF2
                - AfunF1
                - AmacM1
                - AmelC2
                - AmerM2
                - AminM1
                - AquaS1
                - AsinC2
                - AsteI2
                - CpipJ2
                - DmelP6
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

    # open Zarr store
    z = zarr.open(store_location, 'r')
    
    iterable_cut_sites = [(cut_site, z, variation) for cut_site in cut_sites]
    cut_sites = list(tqdm(pool.istarmap(fetch_convar_score, iterable_cut_sites), total=len(iterable_cut_sites), ncols=100))

    return cut_sites