import pandas as pd
import numpy as np
from scipy.stats import rankdata
from sklearn.preprocessing import MinMaxScaler

def rank_guides(df):
    """
        1) Conservation score
        2) MMEJ patterns
            2a) likelihood of MMEJ pattern to arise - MMEJ score
            2b) MMEJ pattern producing OOF deletion
        3) Off-targets
            3a) Mismatch in the seed region & PAM - not likely to be recognized
            3b) Count and score missmatches outside of seed region and PAM
                - MM    0   1   2   3
                - sc    5   3   2   1
        4) Variation - SNPs
            4a) treat SNPs similar to off-targets - in region/out of region - 
                don't effect the score but raise the point

        Weighted sum of all criterias

        Y = x1*w1 + ... xn*wn

        Normalise values on the whole dataset - on all cut_sites


        To maximize target specificity of a sgRNA, we suggest that potential ‘off-target’ 
        genomic sites should be examined by considering the following guidelines: 
            (1) potential ‘off-target’ sites should not be followed by a NGG or NAG PAM sequences; 
            (2) global sequence similarity between sgRNA and potential ‘off-target’ sites should be minimized; 
            (3) single or multiple target mismatch in the “core” sequence is preferred.


    """

    ranking_columns = ['cons_score', 'variants_n', 'sum_score', 'complete_score', 'complete_score', 'offtargets_sum_score']
    neg_ranking_columns = ['sum_score', 'offtargets_sum_score']

    # df.loc[pd.isna(df['offtargets_sum_score']),'offtargets_sum_score'] = 0

    # dfr = df.loc[:,ranking_columns]
    # dfr = dfr.fillna(0.0)
    # dfr.loc[:,neg_ranking_columns] = 1.0 / dfr.loc[:,neg_ranking_columns]
    # dfr.loc[~np.isfinite(dfr[neg_ranking_columns]),neg_ranking_columns] = np.nan

    # print(dfr)

    # scaler = MinMaxScaler()
    # dfr[dfr.columns] = scaler.fit_transform(np.asarray(dfr[dfr.columns], dtype=float))

    # weights = np.array([10, 3, 6, 5, 7])
    # weighted_sum = np.dot(dfr, weights)

    # df['weighted_sum'] = weighted_sum
    # df['rank'] = rankdata(weighted_sum)
    df['weighted_sum'] = np.nan
    df['rank'] = np.nan
    return df