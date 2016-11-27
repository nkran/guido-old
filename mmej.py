import re
import itertools
import math

# TODO
# - look in negative strand
# - microhomology, out-of-frame
# - argv parser - outputs, inputs

# Sequence
seq = 'GATCAGATGCTACGATCGATGTCGCAGCGAGGTGAGGAAGAAAGTGAGGAGGAGGGTGGTAGTGCCACACAGAGAGCTTCGGATG'
max_flanking_length = 40
min_flanking_length = 15
length_weight = 20.0

# Find NGG motives (PAM sites)
# Keep only those which are more than 30 bp downstream m
def find_breaks(sequence):
    breaks_list = []

    pams = [m.start() for m in re.finditer(r'AGG', sequence) if m.start(0) - min_flanking_length > 0 and m.end(0) + min_flanking_length < len(sequence)]
    for pam in pams:

        break_dict = {}

        br = pam - 3
        left = br - max_flanking_length
        right = br + max_flanking_length

        if left < 0:
            left = 0

        # print br, left, right, sequence[pam:pam+3]
        # print br-left, sequence[pam:pam+3], sequence[left:right]
        # print '--'

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
    for k in range(min_kmer_length, len(left_seq)):

        # iterate through sequence
        for i in range(len(left_seq) - k + 1):
            kmer = left_seq[i:i+k]

            # check if k-mer exists in right flanking sequence
            if kmer in right_seq:
                kmers.append(kmer)

    return kmers

# simulates end joining with microhomology
def simulate_end_joining(breaks_list):

    for i, item in enumerate(breaks_list):

        # TODO: remove and enable all breaks
        if i == 4:
            sequence = item['seq']
            br = item['br']

            print sequence

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
                    combination_dict['right'] = '+' * right_deletion_length + right_seq[right_seq_pos:]
                    combination_dict['right_seq'] = right_seq[right_seq_pos:]
                    combination_dict['pattern'] = pattern
                    combination_dict['pattern_score'] = pattern_score
                    combination_dict['deletion_seq'] = deletion_seq
                    combination_dict['frame_shift'] = frame_shift

                    # add to list
                    combination_list.append(combination_dict)

                    # output
                    # print left_seq[:left_seq_pos] + '-' * left_deletion_length + '+' * right_deletion_length + (right_seq[right_seq_pos:]), '\t',
                    # print pattern, '\t', deletion_length,'\t', '\t', frame_shift,'\t', pattern_score,'\t', deletion_seq

            # remove duplicates
            combination_list = [dict(t) for t in set([tuple(sorted(combination_dict.items())) for combination_dict in combination_list])]

            return combination_list

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

    for del_group in sorted(deletion_groups, key=lambda x: x[sorting_field], reverse = True):
        print del_group['score'], '\t', del_group['count'], '\t', del_group['deletion']


# run
break_dict = find_breaks(seq)
pattern_list = simulate_end_joining(break_dict)
group_by_deletion(pattern_list)
