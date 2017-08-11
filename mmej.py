import re
import math
import argparse
import itertools
from Bio.Seq import Seq

parser = argparse.ArgumentParser(description='Microhomology predictor.')
parser.add_argument('--sequence-file', '-i', dest='sequence', help='File with the target sequence.')
parser.add_argument('--good-guides', '-g', dest='good_guides', help='List of good guide sequences.')
parser.add_argument('--bad-guides', '-b', dest='bad_guides', help='List of bad guide sequences.')
parser.add_argument('--max-flanking', '-M', dest='max_flanking_length', help='Max length of flanking region.', default=40)
parser.add_argument('--min-flanking', '-m', dest='min_flanking_length', help='Min length of flanking region.', default=25)
parser.add_argument('--length-weight', '-w', dest='length_weight', help='Length weight - used in scoring.', default=20.0)
parser.add_argument('--rank', help="Output rangking of guides.", default=True)

args = parser.parse_args()

# TODO
# -- remove cut sites with known off target effects
# -- specify custom PAM sequence
# -- allow different outputs
# -- handle FASTA as sequence input

# sequence
max_flanking_length = args.max_flanking_length
min_flanking_length = args.min_flanking_length
length_weight = args.length_weight

with open(args.sequence, 'r') as f:
    seq = f.readline().strip().upper()

with open(args.good_guides, 'r') as f:
    guides_ok = [g.strip() for g in f.readlines()]

with open(args.bad_guides, 'r') as f:
    guides_bad = [g.strip() for g in f.readlines()]


# find NGG motives (PAM sites)
# keep only those which are more than 30 bp downstream m
def find_breaks(sequence):
    breaks_list = []

    pams = [m.start() for m in re.finditer(r'GG', sequence) if m.start(0) - min_flanking_length > 0 and m.end(0) + min_flanking_length < len(sequence)]
    for pam in pams:

        break_dict = {}

        br = pam - 4
        left = br - max_flanking_length
        right = br + max_flanking_length

        if left < 0:
            left = 0

        break_dict['br'] = br - left
        break_dict['seq'] = sequence[left:right]
        break_dict['pam'] = sequence[pam:pam+3]

        breaks_list.append(break_dict)

    return breaks_list


# start with predifined k-mer length and extend it until it finds more
# than one match in sequence
def find_microhomologies(left_seq, right_seq):

    # kmers list
    kmers = []

    # k-mer starting length
    min_kmer_length = 2

    # expand k-mer length
    for k in reversed(xrange(min_kmer_length, len(left_seq))):

        # iterate through sequence
        for i in range(len(left_seq) - k + 1):
            kmer = left_seq[i:i+k]

            if kmer in right_seq:
                kmers.append(kmer)

    return kmers

# simulates end joining with microhomology
def simulate_end_joining(breaks_list):

    cut_sites = []

    for i, item in enumerate(breaks_list):

        sequence = item['seq']
        br = item['br']
        guide = sequence[br-17:br+6]

        # create list for storing combinations
        combination_list = []

        # split sequence at the break
        left_seq = sequence[:br]
        right_seq = sequence[br:]

        # find patterns in both sequences
        patterns = find_microhomologies(left_seq, right_seq)

        # iterate through patterns
        for pattern in patterns:
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
                combination_dict = {}

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
                combination_dict['left'] = left_seq[:left_seq_pos] + '-' * left_deletion_length
                combination_dict['left_seq'] = left_seq[:left_seq_pos]
                combination_dict['left_seq_position'] = left_seq_pos
                combination_dict['right'] = '+' * right_deletion_length + right_seq[right_seq_pos:]
                combination_dict['right_seq'] = right_seq[right_seq_pos:]
                combination_dict['right_seq_position'] = left_seq_pos + len(deletion_seq)
                combination_dict['pattern'] = pattern
                combination_dict['pattern_score'] = pattern_score
                combination_dict['deletion_seq'] = deletion_seq
                combination_dict['frame_shift'] = frame_shift
                combination_dict['guide'] = guide
                combination_dict['sequence'] = sequence

                # add to list
                combination_list.append(combination_dict)

        # remove duplicates and sub microhomologies
        combination_list_filtered = []

        combination_list = [dict(t) for t in set([tuple(sorted(combination_dict.items())) for combination_dict in combination_list])]
        for combination_dict in sorted(combination_list, key=lambda x: x['pattern_score']):

            pass_array = []

            # iterate over previously saved patterns
            for x in combination_list:

                # if the pattern is substring of any previous pattern
                if combination_dict['pattern'] in x['pattern'] and combination_dict['pattern'] != x['pattern']:

                    # get offset position of substring
                    offset_pattern = x['pattern'].find(combination_dict['pattern'])
                    offset_left = combination_dict['left_seq_position'] - (x['left_seq_position'] + offset_pattern)
                    offset_right = combination_dict['right_seq_position'] - (x['right_seq_position'] + offset_pattern)

                    if offset_left == 0 and offset_right == 0:
                        pass_array.append(False)
                    else:
                        pass_array.append(True)
                else:
                    pass_array.append(True)

            # keep only unique mh combinations
            if all(pass_array):
                combination_list_filtered.append(combination_dict)

        cut_sites.append(combination_list_filtered)

    return cut_sites

# groups all combinations of microhomogogies by deletion sequence and sums up the scores for same deletions
def group_by_deletion(pattern_list, sorting_field = 'score'):

    # group list
    deletion_groups = []

    for deletion, group in itertools.groupby(sorted(pattern_list, key=lambda x: x['deletion_seq']), key=lambda x: x['deletion_seq']):

        del_group = {}
        del_group['score'] = 0

        for i, combination in enumerate(group):
            del_group['score'] += combination['pattern_score']

        del_group['count'] = i + 1
        del_group['deletion'] = deletion

        deletion_groups.append(del_group)

    # for del_group in sorted(deletion_groups, key=lambda x: x[sorting_field], reverse = True):
    #     continue
    #     # print del_group['score'], '\t', del_group['count'], '\t', del_group['deletion']

    return deletion_groups

# evaluate guides in specified locus
def evaluate_guides(cut_sites, guides_ok, guides_bad, strand, n_patterns):

    guides = []

    for pattern_list in cut_sites:

        guide = {}

        score = 0
        oof_score = 0

        # sort by pattern score
        sorted_pattern_list = sorted(pattern_list, key=lambda x: x['pattern_score'], reverse=True)[:n_patterns]
        guide_seq = sorted_pattern_list[0]['guide']

        # calculate scores for MH in cut site
        for pattern_dict in sorted_pattern_list:
            if pattern_dict['frame_shift'] == "+":
                oof_score += pattern_dict['pattern_score']

            score += pattern_dict['pattern_score']

        complete_score = oof_score / score * 100

        guide['seq'] = guide_seq
        guide['score'] = complete_score
        guide['strand'] = strand
        guide['patterns'] = sorted_pattern_list

        # check if guide is ok - external checking
        if guide_seq in guides_ok:
            guide['status'] = 'ok'
        elif guide_seq in guides_bad:
            guide['status'] = 'bad'
        else:
            guide['status'] = 'NA'

        guides.append(guide)

    return guides


# run ---------------------

# positive strand
breaks_list_pos = find_breaks(seq)
pattern_list_pos = simulate_end_joining(breaks_list_pos)

# negative strand
revco_seq = Seq(seq).reverse_complement()

breaks_list_neg = find_breaks(str(revco_seq))
pattern_list_neg = simulate_end_joining(breaks_list_neg)

guides_pos = evaluate_guides(pattern_list_pos, guides_ok, guides_bad, "+", 10)
guides_neg = evaluate_guides(pattern_list_neg, guides_ok, guides_bad, "-", 10)

guides_all = guides_neg + guides_pos

if args.rank:
    for guide_dict in sorted(guides_all, key=lambda x: x['score'], reverse=True):

        a = guide_dict['patterns'][0]

        # print a['sequence']

        # print a['left_seq'] + len(a['deletion_seq'])*'-' + a['right_seq']
        # print a['right_seq_position'] * "+"
        # print '\n\n'

        print guide_dict['seq'], '\t', guide_dict['strand'], '\t', guide_dict['status'], '\t', guide_dict['score']
        for pattern in guide_dict['patterns']:

            print pattern['pattern_score'], '\t', len(pattern['deletion_seq']), '\t', pattern['frame_shift'], '\t', pattern['pattern'], '\t', pattern['deletion_seq']

        print "...................................................................."

