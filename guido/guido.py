import os
import re
import sys
import argparse
import pickle

import gffutils
import allel
import multiprocessing as mp
from tqdm import tqdm
from Bio import SeqIO

import guido.log as log
from guido.mmej import simulate_end_joining
from guido.output import render_output
from guido.off_targets import run_bowtie
from guido.convar import apply_conservation_variation_score
from guido.helpers import istarmap, rev_comp, geneset_to_pandas

logger = log.createCustomLogger('root')
ROOT_PATH = os.path.abspath(os.path.dirname(__file__))


def fill_dict(sequence, pams, pam_len, max_flanking_length, region):
    '''
    Creates a list of dictionaries and returns a dataframe of all PAMs with information about:
    position ('break'), pam sequence ('pam'), MMEJ search window ('rel_break' / 'seq'), gRNA ('guide'), and strand ('strand')
    '''

    cuts = []

    for pam_loc in pams:

        cut_dict = {}
        r_chrom, r_start, r_end, r_strand = region
        r_length = r_end - r_start

        # relative cut from the start of the sequence from region obj
        # always 5' -> 3'
        relative_cut_pos = pam_loc - 3
        left_slice = relative_cut_pos - max_flanking_length
        right_slice = relative_cut_pos + max_flanking_length

        if left_slice < 0:
            left_slice = 0
        if right_slice > r_length:
            right_slice = r_length

        seq_slice = slice(left_slice, right_slice)

        if r_strand == '+':
            absolute_cut_pos = r_start + pam_loc - 3
            guide_location = (r_chrom, absolute_cut_pos - 17, absolute_cut_pos + 3 + pam_len)

        if r_strand == '-':
            absolute_cut_pos = r_end - pam_loc + 3
            guide_location = (r_chrom, absolute_cut_pos - 3 - pam_len, absolute_cut_pos + 17)

        # define dict structure -------------------------------
        cut_dict = {
            'pam':                  sequence[pam_loc:pam_loc+pam_len],
            'absolute_cut_pos':     absolute_cut_pos,
            'relative_cut_pos':     relative_cut_pos,
            'relative_cut_pos_seq': relative_cut_pos - left_slice,
            'seq':                  sequence[seq_slice],
            'guide':                sequence[pam_loc-20:pam_loc+pam_len],
            'guide_loc':            guide_location,
            'region':               region,
            'strand':               r_strand,
            'annotation':           [],
            'mmej_patterns':        [],
            'complete_score':       0,
            'sum_score':            0,
            'cons_score':           0,
            'variants':             {},
            'variants_zipped':      {},
            'variants_n':           0,
            'mm':                   {},
            'offtargets_str':       '',
            'offtargets_n':         0,
        }

        if 'N' not in cut_dict['guide']:
            cuts.append(cut_dict)

    return cuts


def find_breaks(region, min_flanking_length, max_flanking_length, pam):
    '''
    Finds Cas9-specific PAM motifs on both strands of a given sequence
    Assumes SpCas9 / 'NGG'-motif by default
    Keeps only those which are more than 30 bp downstream
    '''

    chromosome, start, end, chr_seq = region
    seq = str(chr_seq[start:end].upper())

    iupac_dict = {'A':'A',
                  'C':'C',
                  'G':'G',
                  'T':'T',
                  'R':'[AG]',
                  'Y':'[CT]',
                  'S':'[GC]',
                  'W':'[AT]',
                  'K':'[GT]',
                  'M':'[AC]',
                  'B':'[CGT]',
                  'D':'[AGT]',
                  'H':'[ACT]',
                  'V':'[ACG]',
                  'N':'[ACGT]'}
    iupac_pam = ''.join([iupac_dict[letter] for letter in pam])

    rev_seq = rev_comp(seq)
    pams = [m.start() for m in re.finditer(r'(?=(%s))' % iupac_pam, seq) if m.start(0) - min_flanking_length > 0 and m.end(0) + min_flanking_length < len(seq)]
    rev_pams = [m.start() for m in re.finditer(r'(?=(%s))' % iupac_pam, rev_seq) if m.start(0) - min_flanking_length > 0 and m.end(0) + min_flanking_length < len(rev_seq)]
    pam_len = len(pam)

    cuts_pos = fill_dict(seq, pams, pam_len, max_flanking_length, (chromosome, start, end, '+'))
    cuts_neg = fill_dict(rev_seq, rev_pams, pam_len, max_flanking_length, (chromosome, start, end, '-'))
    cut_sites = cuts_pos + cuts_neg
    cut_sites = sorted(cut_sites, key=lambda x: x['guide_loc'])

    return cut_sites


def parse_args():
    parser = argparse.ArgumentParser(description='Microhomology predictor.')

    parser.add_argument('--sequence-file', '-i', dest='sequence', help='File with the target sequence (TXT or FASTA).')
    parser.add_argument('--region', '-r', dest='region', help='Region in AgamP4 genome [2L:1530-1590].')
    parser.add_argument('--gene', '-G', dest='gene', help='Genome of interest (AgamP4.7 geneset).')
    parser.add_argument('--variants', '-v', dest='variation_store', help='VCF file with variants.', default=False)
    parser.add_argument('--conservation', '-c', dest='conservation_store', help='Path to Zarr store with conservation data.', default=False)
    parser.add_argument('--pam', '-P', dest='pam', help='Protospacer adjacent motif (IUPAC format)', default='NGG')
    parser.add_argument('--threads', '-t', dest='n_threads', type=int, help='Number of threads used.', default=1)
    parser.add_argument('--max-flanking', '-M', type=int, dest='max_flanking_length', help='Max length of flanking region.', default=40)
    parser.add_argument('--min-flanking', '-m', type=int, dest='min_flanking_length', help='Min length of flanking region.', default=25)
    parser.add_argument('--length-weight', '-w', type=float, dest='length_weight', help='Length weight - used in scoring.', default=20.0)
    parser.add_argument('--max-offtargets', type=int, dest='max_offtargets', help='Max number of reported offtargets', default=100)
    parser.add_argument('--n-patterns', '-p', type=int, dest='n_patterns', help='Number of MH patterns used in guide evaluation.', default=5)
    parser.add_argument('--disable-mmej', dest='disable_mmej', type=bool, nargs='?', const=True, help="Disable MMEJ prediction.", default=False)
    parser.add_argument('--disable-off-targets', dest='disable_offtargets', type=bool, nargs='?', const=True, help="Disable off-targets search.", default=False)
    parser.add_argument('--output-folder', '-o', dest='output_folder', help="Output folder.")
    parser.add_argument('--feature-type', '-f', dest='feature', help='Type of genomic feature to focus guide search on.', default=None)
    parser.add_argument('--dump', dest='dump', help="Dump pickled cut_sites object to the output folder.", default=False)

    return parser.parse_args()


def define_genomic_region(chromosome, start, end):
    records = SeqIO.parse(os.path.join(ROOT_PATH, 'data', 'references', 'AgamP4.fa'), "fasta")
    reference = {}

    for record in records:
        reference[record.id] = record

    chr_seq = str(reference[chromosome].seq.upper())
    region = (chromosome, start, end, chr_seq)

    return region


def annotate_guides(cut_site, ann_db, feature):
    '''
    Use GFF annotation to annotate guides with exon name
    (Optional - feature) Only keep guides that are located on a given type of genomic feature
    '''

    location = cut_site['guide_loc']
    pdFeatures = ann_db.query(f'seqid == {repr(location[0])} & start <= {location[2]} & end >= {location[1]}')

    if feature is not None:
        feature_types = [feature for feature in pdFeatures.get('type')]

        if feature == 'intergenic' and 'gene' in feature_types:
            return None

        elif feature == 'intron' and 'gene' not in feature_types or feature == 'intron' and 'exon' in feature_types:
            return None

        else:
            exon_name = ', '.join([exon for exon in set(pdFeatures.get('Name')) if bool(exon)])
            cut_site['annotation'] = exon_name
            return cut_site

    else:
        exon_name = ', '.join([exon for exon in set(pdFeatures.get('Name')) if bool(exon)])
        cut_site['annotation'] = exon_name
        return cut_site

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------
def main():
    ascii_header = r'''

                    ||||||            ||        ||
                  ||    ||                      ||
                  ||        ||    ||  ||    ||||||    ||||
                  ||  ||||  ||    ||  ||  ||    ||  ||    ||
                  ||    ||  ||    ||  ||  ||    ||  ||    ||
                    ||||||    ||||||  ||    ||||||    ||||

                    '''

    print(ascii_header)
    logger.info("Let's dance!")

    args = parse_args()

    max_flanking_length = args.max_flanking_length
    min_flanking_length = args.min_flanking_length
    length_weight = args.length_weight

    # ------------------------------------------------------------
    # Handle input arguments
    # ------------------------------------------------------------

    if args.region or args.gene:
        ann_db = allel.FeatureTable.from_gff3('guido/data/references/AgamP4.7.gff3.gz', attributes=['ID', 'Name'], attributes_fill='')
        ann_db = geneset_to_pandas(ann_db)
    else:
        ann_db = False

    if args.region and args.gene:
        logger.info('Please use only one option for genomic region selection. Use -r or -G.')
        quit()

    if args.region:
        '''
        Option -r: specify the region of interest where guides should be evaluated
        '''

        logger.info('Using AgamP4 reference genome. Region: {}'.format(args.region))

        chromosome = args.region.split(':')[0]
        start = int(args.region.split(':')[1].split('-')[0])
        end = int(args.region.split(':')[1].split('-')[1])

        region = define_genomic_region(chromosome, start, end)

    elif args.gene:
        '''
        Option -G: get genomic region from a gene name
        '''
        logger.info('Using AgamP4 reference genome. Gene: {}'.format(args.gene))

        try:
            gene = ann_db.query(f'ID == {repr(args.gene)}')
        except:
            logger.error('Gene not found: {}'.format(args.gene))
            quit()

        chromosome = ''.join(gene.seqid)
        start = int(gene.get('start'))
        end = int(gene.get('end'))

        region = define_genomic_region(chromosome, start, end)

    elif args.sequence:
        '''
        Option -i: read sequence from FASTA or txt file
        '''
        try:
            record = SeqIO.parse(args.sequence, "fasta").next()
            logger.info('Reading FASTA: {}'.format(args.sequence))
            seq = str(record.seq.upper())
        except:
            logger.info('Reading text sequence', args.sequence)
            with open(args.sequence, 'r') as f:
                seq = f.readline().strip().upper()

        region = ("sequence", 1, len(seq), seq)
    else:
        region = ('seq', 0, 0, False)

    if not args.region and not args.sequence and not args.gene:
        logger.error('Please define the region of interest (-r) or provide the sequence (-i). Use -h for help.')
        quit()

    if args.feature is not None and args.feature != 'intergenic' and args.feature != 'intron':
        feature_types = set([ft for ft in ann_db.get('type')])
        if args.feature not in feature_types:
            logger.error('No feature of this type detected. Features present in current database are the following: intergenic, intron, {}.'.format(', '.join(feature_types)))
            quit()

    if not args.output_folder:
        logger.error('No output folder selected. Please define it by using -o option.')
        quit()
    else:

        # ------------------------------------------------------------
        # Execute
        # ------------------------------------------------------------

        logger.info('Guido is dancing with {} threads ...'.format(args.n_threads))
        pool = mp.Pool(args.n_threads)

        logger.info('Analysing sequence ({} bp) ...'.format(region[2] - region[1]))
        cut_sites = find_breaks(region, min_flanking_length, max_flanking_length, args.pam)

        if ann_db is not False:
            logger.info('Annotating ...')
            iterable_cut_sites = [(cut_site, ann_db, args.feature) for cut_site in cut_sites]
            cut_sites = list(pool.starmap(annotate_guides, iterable_cut_sites))

        if not args.disable_mmej:
            logger.info('Simulating MMEJ ...')
            iterable_cut_sites = [(cut_site, args.n_patterns) for cut_site in cut_sites if cut_site is not None]
            if len(iterable_cut_sites) == 0:
                logger.error('There are no guides that suit your requirements. Try broadening your search parameters.')
                quit()
            else:
                cut_sites = list(pool.starmap(simulate_end_joining, iterable_cut_sites))

        if args.conservation_store or args.variation_store:
            logger.info('Analysing conservation and variation in guides ...')
            cut_sites = apply_conservation_variation_score(cut_sites, args.conservation_store, args.variation_store, pool)

        if not args.disable_offtargets:
            logger.info('Finding offtargets ...')
            cut_sites, targets_df = run_bowtie(cut_sites, args.max_offtargets, args.n_threads)

        pool.close()

        logger.info('Preparing output files ...')
        if not os.path.exists(args.output_folder):
            os.makedirs(args.output_folder)

        if not args.disable_offtargets:
            render_output(cut_sites, args.output_folder, targets_df=targets_df)
        else:
            render_output(cut_sites, args.output_folder)

        if args.dump:
            with open(os.path.join(args.output_folder, 'guides_all.pickle'), 'wb') as f:
                pickle.dump(cut_sites, f, protocol=pickle.HIGHEST_PROTOCOL)
