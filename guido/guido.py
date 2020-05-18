import os
import re
import pickle
import itertools

import pyranges
import multiprocessing as mp
from pyfaidx import Fasta

import guido.log as log
from guido.mmej import simulate_end_joining
from guido.output import render_output
from guido.off_targets import run_bowtie
from guido.convar import apply_conservation_variation_score
from guido.helpers import parse_args, rev_comp, Region

logger = log.createCustomLogger('root')


def fill_dict(sequence, pams, pam_len, max_flanking_length, region, strand):
    '''
    Creates a list of dictionaries and returns a dataframe of all PAMs with
    information about: position ('break'), pam sequence ('pam'), MMEJ search
    window ('rel_break' / 'seq'), gRNA ('guide'), and strand ('strand')
    '''

    cuts = []

    for pam_loc in pams:

        cut_dict = {}
        chrom, start, end, seq, annotation = region
        length = end - start

        # relative cut from the start of the sequence from region obj
        # always 5' -> 3'
        relative_cut_pos = pam_loc - 3
        left_slice = relative_cut_pos - max_flanking_length
        right_slice = relative_cut_pos + max_flanking_length

        if left_slice < 0:
            left_slice = 0
        if right_slice > length:
            right_slice = length

        seq_slice = slice(left_slice, right_slice)

        if strand == '+':
            absolute_cut_pos = start + pam_loc - 3
            guide_location = (
                chrom,
                absolute_cut_pos - 17,
                absolute_cut_pos + 3 + pam_len,
            )

        if strand == '-':
            absolute_cut_pos = end - pam_loc + 3
            guide_location = (
                chrom,
                absolute_cut_pos - 3 - pam_len,
                absolute_cut_pos + 17,
            )

        guide_annotations = []
        if annotation is not None:
            guide_annotations = annotation.query(f'((Start <= {absolute_cut_pos}) & (End >= {absolute_cut_pos}))')

        # define dict structure -------------------------------
        cut_dict = {
            'pam': sequence[pam_loc:pam_loc+pam_len],
            'absolute_cut_pos': absolute_cut_pos,
            'relative_cut_pos': relative_cut_pos,
            'relative_cut_pos_seq': relative_cut_pos - left_slice,
            'seq': sequence[seq_slice],
            'guide': sequence[pam_loc-20:pam_loc+pam_len],
            'guide_loc': guide_location,
            'strand': strand,
            'annotation': guide_annotations,
            'mmej_patterns': [],
            'complete_score': 0,
            'sum_score': 0,
            'cons_score': 0,
            'variants': {},
            'variants_zipped': {},
            'variants_n': 0,
            'mm': {},
            'offtargets_str': '',
            'offtargets_n': 0,
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

    chromosome, start, end, seq, annotation = region
    iupac_dict = {
        'A': 'A',
        'C': 'C',
        'G': 'G',
        'T': 'T',
        'R': '[AG]',
        'Y': '[CT]',
        'S': '[GC]',
        'W': '[AT]',
        'K': '[GT]',
        'M': '[AC]',
        'B': '[CGT]',
        'D': '[AGT]',
        'H': '[ACT]',
        'V': '[ACG]',
        'N': '[ACGT]',
    }
    iupac_pam = ''.join([iupac_dict[letter] for letter in pam])

    rev_seq = rev_comp(seq)
    pams = [
        m.start()
        for m in re.finditer(r'(?=(%s))' % iupac_pam, seq)
        if m.start(0) - min_flanking_length > 0
        and m.end(0) + min_flanking_length < len(seq)
    ]
    rev_pams = [
        m.start()
        for m in re.finditer(r'(?=(%s))' % iupac_pam, rev_seq)
        if m.start(0) - min_flanking_length > 0
        and m.end(0) + min_flanking_length < len(rev_seq)
    ]
    pam_len = len(pam)

    cuts_pos = fill_dict(seq, pams, pam_len, max_flanking_length, region, '+')
    cuts_neg = fill_dict(rev_seq, rev_pams, pam_len, max_flanking_length, region, '-')
    cut_sites = cuts_pos + cuts_neg
    cut_sites = sorted(cut_sites, key=lambda x: x['guide_loc'])

    return cut_sites


def define_genomic_region(chromosome, start, end, genome_file, annotation_db=None):
    genome = Fasta(genome_file)
    region_seq = genome[chromosome][start:end].seq.upper()

    if annotation_db is not None:
        annotation = annotation_db.query(
            f'(Chromosome == {repr(chromosome)}) &  \
              (((Start >= {start - 1}) & (Start <= {end + 1})) | \
              ((End >= {start - 1}) & (End <= {end + 1})))')
    else:
        annotation = None

    return Region(
        chromosome=chromosome,
        start=start,
        end=end,
        sequence=region_seq,
        annotation=annotation,
    )


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

    # ------------------------------------------------------------
    # Handle input arguments
    # ------------------------------------------------------------
    if not args.region and not args.sequence and not args.gene:
        logger.error(
            'Please define the region of interest (-r) or provide the sequence (-i). Use -h for help.'
        )
        quit()

    if args.region and args.gene:
        logger.info(
            'Please use only one option for genomic region selection. Use -r or -G.'
        )
        quit()

    if not args.genome_info:
        logger.error(
            "Please provide the genome information file destination using -g argument."
        )
    else:
        genome_info = pickle.load(open(args.genome_info, "rb"))
        genome_index_path = genome_info['genome_index_path']
        genome_file_path = genome_info['genome_file']
        annotation_file_path = genome_info['annotation_file']
        ann_ext = genome_info['ann_ext']
        if annotation_file_path is None:
            if args.gene or args.feature:
                logger.error('Please provide an Annotation file while using guido-build if you wish to use Gene and/or Feature arguments.')
                quit()

    if args.region or args.gene:
        if annotation_file_path is not None:
            if ann_ext in ['.gff3']:
                ann_db = pyranges.read_gff3(genome_info['annotation_file'], as_df=True)
                ann_db.rename(columns={'Name':'Exon'}, inplace=True) # TODO - Keep only number after the final 'E' as exon_number
            elif ann_ext in ['.gtf']:
                ann_db = pyranges.read_gtf(genome_info['annotation_file'], as_df=True)
                ann_db.rename(columns={'gene_id':'ID','exon_number':'Exon'}, inplace=True)
        else:
            ann_db = None

    if args.feature and len(ann_db) > 0:
        feature_types = ann_db['Feature'].unique()
        if args.feature not in feature_types:
            logger.error(
                f'No feature of this type detected. Features present in current database are the following: {list(feature_types)}'
            )
            quit()

    if args.region:
        '''
        Option -r: specify the region of interest where guides should be evaluated
        '''

        logger.info('Using AgamP4 reference genome. Region: {}'.format(args.region))

        chromosome = args.region.split(':')[0]
        start = int(args.region.split(':')[1].split('-')[0])
        end = int(args.region.split(':')[1].split('-')[1])

    elif args.gene:
        '''
        Option -G: get genomic region from a gene name
        '''
        logger.info('Using AgamP4 reference genome. Gene: {}'.format(args.gene))

        try:
            gene = ann_db.query(f'ID == {repr(args.gene)} & Feature == "gene"')
        except Exception:
            logger.error('Gene not found: {}'.format(args.gene))
            quit()

        chromosome = gene.Chromosome.values[0]
        start = int(gene.Start)
        end = int(gene.End)

    if args.feature:
        overlapping_features = ann_db.query(
            f'(Feature == {repr(args.feature)}) & \
                                                (Chromosome == {repr(chromosome)}) &  \
                                                (((Start >= {start}) & (Start <= {end})) | \
                                                ((End >= {start}) & (End <= {end})))'
        )
        if len(overlapping_features) > 0:
            regions = [
                define_genomic_region(chromosome, start, end, genome_file_path, ann_db)
                for chromosome, start, end in overlapping_features[
                    ['Chromosome', 'Start', 'End']
                ].values
            ]
        else:
            logger.error('No feature of this type detected in the provided region.')
            quit()
    else:
        regions = [define_genomic_region(chromosome, start, end, genome_file_path, ann_db)]

    if not args.output_folder:
        logger.error('No output folder selected. Please define it by using -o option.')
        quit()
    else:

        # ------------------------------------------------------------
        # Execute
        # ------------------------------------------------------------

        logger.info('Guido is dancing with {} threads ...'.format(args.n_threads))
        pool = mp.Pool(args.n_threads)

        logger.info('Analysing sequence ...')

        cut_sites = [
            find_breaks(region, min_flanking_length, max_flanking_length, args.pam)
            for region in regions
        ]
        cut_sites = list(itertools.chain(*cut_sites))

        if not args.disable_mmej:
            logger.info('Simulating MMEJ ...')
            iterable_cut_sites = [(cut_site, args.n_patterns) for cut_site in cut_sites]
            if not iterable_cut_sites:
                logger.error(
                    'There are no guides that suit your requirements. Try broadening your search parameters.'
                )
                quit()
            else:
                cut_sites = pool.starmap(simulate_end_joining, iterable_cut_sites)

        if args.conservation_store or args.variation_store:
            logger.info('Analysing conservation and variation in guides ...')
            cut_sites = apply_conservation_variation_score(
                cut_sites, args.conservation_store, args.variation_store, pool
            )

        if not args.disable_offtargets:
            logger.info('Finding offtargets ...')
            cut_sites, targets_df = run_bowtie(
                cut_sites, args.max_offtargets, genome_index_path, args.n_threads
            )

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
