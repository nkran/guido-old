import re
import math
import itertools
import pandas as pd
import guido.log as log

logger = log.createCustomLogger('MMEJ')


def find_microhomologies(left_seq, right_seq):
    '''
    Start with predefined k-mer length and extend it until it finds more
    than one match in sequence.
    '''

    # kmers list
    kmers = []

    # k-mer starting length
    min_kmer_length = 2

    # expand k-mer length
    for k in reversed(range(min_kmer_length, len(left_seq))):

        # iterate through sequence
        for i in range(len(left_seq) - k + 1):
            kmer = left_seq[i:i+k]

            if kmer in right_seq and kmer not in kmers:
                kmers.append(kmer)

    return kmers


def generate_mmej_patterns(rel_break_position, sequence, length_weight):
    '''
    Generates MMEJ patterns for the cut site and returns a dictionary of patterns
    '''

    # create list for storing MH patterns
    mmej_patterns = []

    if 'N' not in sequence:
        # split sequence at the break
        left_seq = sequence[:rel_break_position]
        right_seq = sequence[rel_break_position:]

        # find patterns in both sequences
        sequence_patterns = find_microhomologies(left_seq, right_seq)

        # iterate through patterns
        for pattern in sequence_patterns:
            p = re.compile(pattern)

            # find positions of patterns in each sequence
            left_positions = [m.start() for m in p.finditer(left_seq)]
            right_positions = [m.start() for m in p.finditer(right_seq)]

            # GC count for pattern
            pattern_GC = len(re.findall('G', pattern)) + len(re.findall('C', pattern))

            # get combinations
            pos_combinations = list(itertools.product(left_positions, right_positions))

            # generate microhomology for every combination
            for combination in pos_combinations:
                # save output to dict
                pattern_dict = {}

                # left side
                left_seq_pos = combination[0]
                left_deletion_length = len(left_seq) - left_seq_pos

                # right side
                right_seq_pos = combination[1]
                right_deletion_length = right_seq_pos

                # deletion length and sequence
                deletion_length = left_deletion_length + right_deletion_length
                deletion_seq = left_seq[left_seq_pos:] + right_seq[:right_seq_pos]

                # score pattern
                length_factor =  round(1 / math.exp((deletion_length) / (length_weight)), 3)
                pattern_score = 100 * length_factor * ((len(pattern) - pattern_GC) + (pattern_GC * 2))

                # frame shift
                if deletion_length % 3 == 0:
                    frame_shift = "-"
                else:
                    frame_shift = "+"

                # create dictionary
                pattern_dict['left'] = left_seq[:left_seq_pos] + '-' * left_deletion_length
                pattern_dict['left_seq'] = left_seq[:left_seq_pos]
                pattern_dict['left_seq_position'] = left_seq_pos
                pattern_dict['right'] = '+' * right_deletion_length + right_seq[right_seq_pos:]
                pattern_dict['right_seq'] = right_seq[right_seq_pos:]
                pattern_dict['right_seq_position'] = left_seq_pos + len(deletion_seq)
                pattern_dict['pattern'] = pattern
                pattern_dict['pattern_len'] = len(pattern)
                pattern_dict['pattern_score'] = pattern_score
                pattern_dict['deletion_seq'] = deletion_seq
                pattern_dict['frame_shift'] = frame_shift

                # add to list
                mmej_patterns.append(pattern_dict)

    return mmej_patterns


def simulate_end_joining(cut_site, n_patterns):
    '''
    Simulates end joining with microhomologies
    '''
    
    length_weight = 20
    column_dtypes = {
        'break': 'uint64',
        'break_abs': 'uint64',
        'rel_break': 'uint16',
        'seq': str,
        'pam': str,
        'guide': str,
        'strand': str,
        'mmej_patterns': 'object',
    }

    mmej_patterns = generate_mmej_patterns(cut_site['rel_break'], cut_site['seq'], length_weight)
    
    if mmej_patterns:
        sorted_mmej_patterns = pd.DataFrame(mmej_patterns).drop_duplicates() \
                                                          .sort_values('pattern_score') \
                                                          .groupby(['deletion_seq', 'left_seq_position', 'right_seq_position']) \
                                                          .first() \
                                                          .reset_index() \
                                                          .sort_values('pattern_score', ascending=False)

        sorted_mmej_patterns = sorted_mmej_patterns.head(n_patterns)

        oof_score = sorted_mmej_patterns[sorted_mmej_patterns['frame_shift'] == '+'].loc[:,'pattern_score'].sum()
        score = sorted_mmej_patterns.loc[:,'pattern_score'].sum()
        complete_score = oof_score / score * 100

        cut_site['mmej_patterns'] = sorted_mmej_patterns
        cut_site['complete_score'] = complete_score
        cut_site['sum_score'] = score
        cut_site['mmej_patterns'] = sorted_mmej_patterns
    else:
        cut_site['mmej_patterns'] = []
        cut_site['complete_score'] = 0
        cut_site['sum_score'] = 0
        cut_site['mmej_patterns'] = 0
    
    return cut_site