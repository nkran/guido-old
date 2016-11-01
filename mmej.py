import re
import Bio
import itertools

# TODO
# - look in negative strand
# - score patterns, microhomology, out-of-frame

#sequence
seq = 'CGAGCGCGGGGCAGGTGCCCGCTGGAACTCGCGCCTCGCAGCGCTGGGCGGCCGGGGCCGGGCAGGGTAGTGCGGGAAGATCGGGGTCTGGGGTCGGTGCCGGCGGGACTCCGAAAGGAGGGAGCCGGG'
max_flanking_length = 40
min_flanking_length = 15


# Find NGG motives (PAM sites)
# Keep only those which are more than 30 bp downstream
def find_breaks(sequence):
    breaks_list = []

    pams = [m.start() for m in re.finditer(r'[ACGT]GG', sequence) if m.start(0) - min_flanking_length > 0 and m.end(0) + min_flanking_length < len(sequence)]
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


def simulate_end_joining(breaks_list):

    for i, item in enumerate(breaks_list):
        if i == 0:
            sequence = item['seq']
            br = item['br']

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

                # get combinations
                pos_combinations = list(itertools.product(left_positions, right_positions))

                for combination in pos_combinations:

                    left_seq_pos = combination[0]
                    left_deletion_length = len(left_seq) - left_seq_pos

                    right_seq_pos = combination[1]
                    right_deletion_length = right_seq_pos

                    print 'microhomology:', pattern
                    print 'wild type:\n', sequence
                    print left_seq[:left_seq_pos] + '-' * left_deletion_length + '+' * right_deletion_length + (right_seq[right_seq_pos:])
                    print '\n'


break_dict = find_breaks(seq)
simulate_end_joining(break_dict)

