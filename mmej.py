from __future__ import print_function
import os
import re
import math
import argparse
import itertools
import vcf
from Bio.Seq import Seq
from Bio import SeqIO


parser = argparse.ArgumentParser(description='Microhomology predictor.')
parser.add_argument('--sequence-file', '-i', dest='sequence', help='File with the target sequence (TXT or FASTA).')
parser.add_argument('--region', '-r', dest='region', help='Region in AgamP4 genome [2L:1530-1590].')
parser.add_argument('--variants', '-v', dest='variants', help='VCF file with variants.')
parser.add_argument('--good-guides', '-g', dest='good_guides', help='List of good guide sequences.')
parser.add_argument('--bad-guides', '-b', dest='bad_guides', help='List of bad guide sequences.')
parser.add_argument('--max-flanking', '-M', dest='max_flanking_length', help='Max length of flanking region.', default=40)
parser.add_argument('--min-flanking', '-m', dest='min_flanking_length', help='Min length of flanking region.', default=25)
parser.add_argument('--length-weight', '-w', dest='length_weight', help='Length weight - used in scoring.', default=20.0)
parser.add_argument('--n-patterns', '-p', dest='n_patterns', help='Number of MH patterns used in guide evaluation.', default=5)
parser.add_argument('--output-folder', '-o', dest='output_folder', help="Output folder.")

args = parser.parse_args()


# find NGG motives (PAM sites)
# keep only those which are more than 30 bp downstream m
def find_breaks(sequence, region):
    breaks_list = []

    pams = [m.start() for m in re.finditer(r'(?=([ACTG]GG))', sequence) if m.start(0) - min_flanking_length > 0 and m.end(0) + min_flanking_length < len(sequence)]

    for pam in pams:

        break_dict = {}

        br = pam - 3
        left = br - max_flanking_length
        right = br + max_flanking_length

        if left < 0:
            left = 0

        break_dict['br'] = br - left
        break_dict['seq'] = sequence[left:right]
        break_dict['pam'] = sequence[pam:pam+3]

        break_dict['left_flank'] = sequence[left]

        break_dict['br_abs'] = region[1] + br
        break_dict['chr'] = region[0]
        break_dict['ref_start'] = region[1]
        break_dict['ref_end'] = region[2]

        breaks_list.append(break_dict)

    return breaks_list

def extract_flanks(chr_seq, break_dict_list):

    break_dicts = []

    for break_dict in break_dict_list:

        br = break_dict['br_abs']

        if br - 2000 > 0:
            left_flank_start = br - 2000
        else:
            left_flank_start = 0

        if br + 2000 < len(chr_seq):
            right_flank_end = br + 2000
        else:
            right_flank_end = len(chr_seq)

        break_dict['left_flank_seq'] = chr_seq[left_flank_start:br]
        break_dict['right_flank_seq'] = chr_seq[br:right_flank_end]

        break_dicts.append(break_dict)

    return break_dicts

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

        cut_site = {}

        br = item['br']
        seq = item['seq']

        cut_site['sequence'] = item['seq']
        cut_site['br'] = item['br']
        cut_site['guide'] = seq[br-17:br+6]
        cut_site['guide_loc'] = item['chr'] + ':' + str(item['br_abs'] - 17) + '-' + str(item['br_abs'] + 6)
        cut_site['left_flank_seq'] = item['left_flank_seq']
        cut_site['right_flank_seq'] = item['right_flank_seq']

        # create list for storing combinations
        cut_site['combination_list'] = []

        # split sequence at the break
        left_seq = seq[:br]
        right_seq = seq[br:]

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

                # add to list
                cut_site['combination_list'].append(combination_dict)

        # remove duplicates and sub microhomologies
        combination_list_filtered = []

        combination_list = [dict(t) for t in set([tuple(sorted(combination_dict.items())) for combination_dict in cut_site['combination_list']])]
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

        cut_site['combination_list'] = combination_list_filtered

        cut_sites.append(cut_site)

    return cut_sites

# groups all combinations of microhomogogies by deletion sequence and sums up the scores for same deletions
def group_by_deletion(pattern_list, sorting_field = 'score'):

    # group list
    deletion_groups = []

    for deletion, group in itertools.groupby(sorted(pattern_list, key=lambda x: x['deletion_seq']), key=lambda x: x['deletion_seq']):

        del_group = {}
        del_group['score'] = 0
        del_group['pattern'] = []
        del_group['guide'] = set()
        for i, combination in enumerate(group):
            del_group['pattern'].append(combination['pattern'])
            del_group['score'] += combination['pattern_score']
            del_group['guide'].add(combination['guide'])

        del_group['count'] = i + 1
        del_group['deletion'] = deletion
        del_group['deletion_len'] = len(deletion)
        deletion_groups.append(del_group)

    for del_group in sorted(deletion_groups, key=lambda x: del_group['deletion_len']):
        print(del_group['score'], '\t', del_group['count'], '\t', del_group['deletion_len'], '\t', del_group['deletion'], '\t', del_group['pattern'], '\t', del_group['guide'])

    return deletion_groups

# evaluate guides in specified locus
def evaluate_guides(cut_sites, guides_ok, guides_bad, strand, n_patterns, var_positions):

    guides = []

    for cut_site in cut_sites:
        guide = {}
        score = 0
        oof_score = 0

        # sort by pattern score
        sorted_pattern_list = sorted(cut_site['combination_list'], key=lambda x: x['pattern_score'], reverse=True)[:n_patterns]
        guide_seq = cut_site['guide']
        guide_loc = cut_site['guide_loc']

        guide_start = int(guide_loc.split(':')[1].split('-')[0])
        guide_end = int(guide_loc.split(':')[1].split('-')[1])

        # calculate SNP penalty
        snp_score = 0

        if var_positions:
            for pos in var_positions:
                if guide_start <= pos and guide_end >= pos:
                    snp_score += 1
        else:
            print("No variants")

        # calculate scores for MH in cut site
        for pattern_dict in sorted_pattern_list:
            if pattern_dict['frame_shift'] == "+":
                oof_score += pattern_dict['pattern_score']

            score += pattern_dict['pattern_score']

        complete_score = oof_score / score * 100

        guide['seq'] = guide_seq
        guide['guide_loc'] = guide_loc
        guide['score'] = complete_score
        guide['snp_score'] = snp_score
        guide['sum_score'] = score
        guide['strand'] = strand
        guide['patterns'] = sorted_pattern_list
        guide['left_flank_seq'] = cut_site['left_flank_seq']
        guide['right_flank_seq'] = cut_site['right_flank_seq']

        # check if guide is ok - external checking
        if guide_seq in guides_ok:
            guide['status'] = 'ok'
        elif guide_seq in guides_bad:
            guide['status'] = 'bad'
        else:
            guide['status'] = 'NA'

        guides.append(guide)

    return guides


# TODO
# -- remove cut sites with known off target effects
# -- specify custom PAM sequence
# -- allow different outputs

# sequence
max_flanking_length = int(args.max_flanking_length)
min_flanking_length = int(args.min_flanking_length)
length_weight = args.length_weight

# sequence input -------------------------------------------------------
if args.sequence:
    try:
        record = SeqIO.parse(args.sequence, "fasta").next()
        print('Reading FASTA', args.sequence)
        seq = str(record.seq.upper())
    except:
        print('Reading text sequence', args.sequence)
        with open(args.sequence, 'r') as f:
            seq = f.readline().strip().upper()

# sequence input from reference genome ---------------------------------
if args.region:
    print('Using AgamP4 sequence')

    records = SeqIO.parse("references/AgamP4.fa", "fasta")
    reference = {}

    chromosome = args.region.split(':')[0]
    start = int(args.region.split(':')[1].split('-')[0])
    end = int(args.region.split(':')[1].split('-')[1])

    for record in records:
        reference[record.id] = record

    seq = str(reference[chromosome].seq[start:end].upper())
    chr_seq = str(reference[chromosome].seq.upper())

    region = (chromosome, start, end)
else:
    region = ('seq', 0, 0)

# good guides ----------------------------------------------------------
if args.good_guides:
    with open(args.good_guides, 'r') as f:
        guides_ok = [g.strip() for g in f.readlines()]
else:
    guides_ok = []

# bad guides -----------------------------------------------------------
if args.bad_guides:
    with open(args.bad_guides, 'r') as f:
        guides_bad = [g.strip() for g in f.readlines()]
else:
    guides_bad = []

# variants input -------------------------------------------------------
if args.variants and args.region:
    with open(args.variants, 'r') as vcf_file:
        vcf_reader = vcf.Reader(vcf_file)

        print('Reading VCF for region {}:{}-{}'.format(chromosome, start, end))
        var_positions = [v.POS for v in vcf_reader.fetch(chromosome, start, end)]


# ----------------------------------------------------------------------
# run
# ----------------------------------------------------------------------

# neg strand
revco_seq = str(Seq(seq).reverse_complement())
revco_chr_seq = str(Seq(chr_seq).reverse_complement())

# positive strand
breaks_list_pos = find_breaks(seq, region)
breaks_list_neg = find_breaks(revco_seq, region)

if args.region:
    breaks_list_pos = extract_flanks(chr_seq, breaks_list_pos)
    breaks_list_neg = extract_flanks(chr_seq, breaks_list_neg)

pattern_list_pos = simulate_end_joining(breaks_list_pos)
pattern_list_neg = simulate_end_joining(breaks_list_neg)

# evaluate guides
guides_pos = evaluate_guides(pattern_list_pos, guides_ok, guides_bad, "+", int(args.n_patterns), var_positions)
guides_neg = evaluate_guides(pattern_list_neg, guides_ok, guides_bad, "-", int(args.n_patterns), var_positions)

# merge guides from positive and negative strand
guides_all = guides_neg + guides_pos

# output
if args.output_folder:

    # create output dir if it doesn't exist
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    # print the list of guides
    with open(args.output_folder + '/guides_list_' + str(args.n_patterns) + '.txt', 'w') as f:

        print("guide_sequence\tgenomic_location\tstrand\toff_target_analysis\tMMEJ_score\tMMEJ_sum_score\tMMEJ_top_score\tMMEJ_out_of_frame_del\tguide_snp_count\tleft_flank\tright_flank", file = f)

        for guide_dict in sorted(guides_all, key=lambda x: (x['score'], x['sum_score']), reverse=True):
            mmej_frames = ('\t'.join('{}'.format(pattern['frame_shift']) for pattern in guide_dict['patterns']))

            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(guide_dict['seq'], guide_dict['guide_loc'], guide_dict['strand'], guide_dict['status'], guide_dict['score'], guide_dict['sum_score'], guide_dict['patterns'][0]['pattern_score'], mmej_frames, guide_dict['snp_score'], guide_dict['left_flank_seq'], guide_dict['right_flank_seq']), file = f)

    # print the list of guides and their MH patterns
    with open(args.output_folder + '/guides_list_mh_' + str(args.n_patterns) + '.txt', 'w') as f:

        for guide_dict in sorted(guides_all, key=lambda x: (x['score'], x['sum_score']), reverse=True):
            print("{}\t{}\t{}\t{}\t{}\t{}".format(guide_dict['seq'], guide_dict['guide_loc'], guide_dict['strand'], guide_dict['status'], guide_dict['score'], guide_dict['sum_score']), file = f)

            for pattern in guide_dict['patterns']:
                print("{}\t{}\t{}\t{}\t{}".format(pattern['pattern_score'], len(pattern['deletion_seq']), pattern['frame_shift'], pattern['pattern'], pattern['deletion_seq']), file = f)

            print("....................................................................", file = f)

