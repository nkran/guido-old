from __future__ import print_function
import os
import re
import math
import argparse
import itertools
import vcf
from Bio.Seq import Seq, reverse_complement
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


def find_breaks(sequence):
    '''
    Find NGG motives (PAM sites)
    Keep only those which are more than 30 bp downstream
    '''

    pams = [m.start() for m in re.finditer(r'(?=([ACTG]GG))', sequence) if m.start(0) - args.min_flanking_length > 0 and m.end(0) + args.min_flanking_length < len(sequence)]
    breaks_list = []

    for pam in pams:

        break_dict = {}

        br = pam - 3
        left = br - args.max_flanking_length
        right = br + args.max_flanking_length

        if left < 0:
            left = 0

        break_dict['br'] = br - left
        break_dict['break'] = br
        break_dict['seq'] = sequence[left:right]
        break_dict['pam'] = sequence[pam:pam+3]
        break_dict['guide'] = sequence[br-17:br+6]

        breaks_list.append(break_dict)

    return breaks_list


def get_cut_sites(region):
    '''
    Get cutsites for positive and negative strand in the context
    of provided genomic region.
    '''

    chromosome, start, end, chr_seq = region

    seq = str(chr_seq[start:end].upper())
    rev_seq = reverse_complement(seq)

    cuts_plus = find_breaks(seq)
    cuts_minus = find_breaks(rev_seq)

    for cut in cuts_plus:
        break_abs = cut['break'] + start

        cut.update({'strand': '+'})
        cut.update({'break_abs': break_abs})
        cut.update({'guide_loc': (chromosome, break_abs - 17, break_abs + 6)})

        cut.update({'left_flank_seq': chr_seq[break_abs - 2000:break_abs]})
        cut.update({'right_flank_seq': chr_seq[break_abs:break_abs + 2000]})

    for cut in cuts_minus:
        break_abs = end - cut['break']

        cut.update({'strand': '-'})
        cut.update({'break_abs': break_abs + 1})
        cut.update({'guide_loc': (chromosome, break_abs - 5, break_abs + 17)})

        cut.update({'left_flank_seq': chr_seq[break_abs - 2000:break_abs]})
        cut.update({'right_flank_seq': chr_seq[break_abs:break_abs + 2000]})

    # merge lists
    cut_list = cuts_plus + cuts_minus

    return cut_list


def find_microhomologies(left_seq, right_seq):
    '''
    Start with predifined k-mer length and extend it until it finds more
    than one match in sequence
    '''
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

def simulate_end_joining(cut_list):
    '''
    Simulates end joining with microhomologies
    '''

    for i, cut in enumerate(cut_list):

        br = cut['br']
        seq = cut['seq']

        # create list for storing MH patterns
        cut.update({'pattern_list': []})

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
                length_factor =  round(1 / math.exp((deletion_length) / (args.length_weight)), 3)
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
                pattern_dict['pattern_score'] = pattern_score
                pattern_dict['deletion_seq'] = deletion_seq
                pattern_dict['frame_shift'] = frame_shift

                # add to list
                cut['pattern_list'].append(pattern_dict)

        # remove duplicates and sub microhomologies
        pattern_list_filtered = []

        pattern_list = [dict(t) for t in set([tuple(sorted(pattern_dict.items())) for pattern_dict in cut['pattern_list']])]
        for pattern_dict in sorted(pattern_list, key=lambda x: x['pattern_score']):

            pass_array = []

            # iterate over previously saved patterns
            for x in pattern_list:

                # if the pattern is substring of any previous pattern
                if pattern_dict['pattern'] in x['pattern'] and pattern_dict['pattern'] != x['pattern']:

                    # get offset position of substring
                    offset_pattern = x['pattern'].find(pattern_dict['pattern'])
                    offset_left = pattern_dict['left_seq_position'] - (x['left_seq_position'] + offset_pattern)
                    offset_right = pattern_dict['right_seq_position'] - (x['right_seq_position'] + offset_pattern)

                    if offset_left == 0 and offset_right == 0:
                        pass_array.append(False)
                    else:
                        pass_array.append(True)
                else:
                    pass_array.append(True)

            # keep only unique mh patterns
            if all(pass_array):
                pattern_list_filtered.append(pattern_dict)

        cut.update({'pattern_list': pattern_list_filtered})

    return cut_sites


def evaluate_guides(cut_sites, guides_ok, guides_bad, n_patterns, var_positions):
    '''
    Score guides and include information about SNPs and out-of-frame deletions
    '''
    guides = []

    for cut_site in cut_sites:
        guide = {}
        score = 0
        oof_score = 0

        # sort by pattern score
        sorted_pattern_list = sorted(cut_site['pattern_list'], key=lambda x: x['pattern_score'], reverse=True)[:n_patterns]
        guide_seq = cut_site['guide']

        chromosome, guide_start, guide_end = cut_site['guide_loc']

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

        cut_site.update({'complete_score': complete_score})
        cut_site.update({'sum_score': score})
        cut_site.update({'snp_score': snp_score})
        cut_site.update({'top_patterns': sorted_pattern_list})


        # check if guide is ok - external checking
        if guide_seq in guides_ok:
            cut_site.update({'status': 'ok'})
        elif guide_seq in guides_bad:
            cut_site.update({'status': 'bad'})
        else:
            cut_site.update({'status': 'NA'})

    return cut_sites


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

if args.region:
    print('Using AgamP4 sequence')

    records = SeqIO.parse("references/AgamP4.fa", "fasta")
    reference = {}

    chromosome = args.region.split(':')[0]
    start = int(args.region.split(':')[1].split('-')[0])
    end = int(args.region.split(':')[1].split('-')[1])

    for record in records:
        reference[record.id] = record

    chr_seq = str(reference[chromosome].seq.upper())
    region = (chromosome, start, end, chr_seq)
else:
    region = ('seq', 0, 0, False)

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


cut_sites = get_cut_sites(region)
cut_sites = simulate_end_joining(cut_sites)
cut_sites = evaluate_guides(cut_sites, guides_ok, guides_bad, int(args.n_patterns), var_positions)

# output
if args.output_folder:

    # create output dir if it doesn't exist
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    # print the list of guides
    with open(args.output_folder + '/guides_list_' + str(args.n_patterns) + '.txt', 'w') as f:

        print("guide_sequence\tgenomic_location\tstrand\toff_target_analysis\tMMEJ_score\tMMEJ_sum_score\tMMEJ_top_score\tMMEJ_out_of_frame_del\tguide_snp_count\tleft_flank\tright_flank", file = f)

        for guide_dict in sorted(cut_sites, key=lambda x: (x['complete_score'], x['sum_score']), reverse=True):
            mmej_frames = ('\t'.join('{}'.format(pattern['frame_shift']) for pattern in guide_dict['top_patterns']))
            location = guide_dict['guide_loc'][0] + ':' + str(guide_dict['guide_loc'][1]) + '-' + str(guide_dict['guide_loc'][2])

            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(guide_dict['guide'], location, guide_dict['strand'], guide_dict['status'], guide_dict['complete_score'], guide_dict['sum_score'], guide_dict['top_patterns'][0]['pattern_score'], mmej_frames, guide_dict['snp_score'], guide_dict['left_flank_seq'], guide_dict['right_flank_seq']), file = f)

    # print the list of guides and their MH patterns
    with open(args.output_folder + '/guides_list_mh_' + str(args.n_patterns) + '.txt', 'w') as f:

        for guide_dict in sorted(cut_sites, key=lambda x: (x['complete_score'], x['sum_score']), reverse=True):
            location = guide_dict['guide_loc'][0] + ':' + str(guide_dict['guide_loc'][1]) + '-' + str(guide_dict['guide_loc'][2])

            print("{}\t{}\t{}\t{}\t{}\t{}".format(guide_dict['guide'], location, guide_dict['strand'], guide_dict['status'], guide_dict['complete_score'], guide_dict['sum_score']), file = f)
            for pattern in guide_dict['top_patterns']:
                print("{}\t{}\t{}\t{}\t{}".format(pattern['pattern_score'], len(pattern['deletion_seq']), pattern['frame_shift'], pattern['pattern'], pattern['deletion_seq']), file = f)

            print("....................................................................", file = f)

