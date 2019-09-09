import os
import re
import sys
import argparse

import vcf
import gffutils
import multiprocessing as mp
from tqdm import tqdm
from Bio import SeqIO

import guido.log as log
from guido.mmej import simulate_end_joining
from guido.off_targets import run_bowtie
from guido.helpers import istarmap, rev_comp


logger = log.createCustomLogger('root')
ROOT_PATH = os.path.abspath(os.path.dirname(__file__))


def fill_dict(sequence, start, pams, pam_len, max_flanking_length, region):
    '''
    Creates a list of dictionaries and returns a dataframe of all PAMs with information about:
    position ('break'), pam sequence ('pam'), MMEJ search window ('rel_break' / 'seq'), gRNA ('guide'), and strand ('strand')
    '''

    cuts = []

    for pam in pams:
        
        cut_dict = {}
        strand = region[3]
        br = pam - pam_len
        left = br - max_flanking_length
        if left < 0:
            left = 0
        right = br + max_flanking_length

        cut_dict['break'] = br
        cut_dict['abs_break'] = br + start
        cut_dict['rel_break'] = br - left
        cut_dict['seq'] = sequence[left:right]
        cut_dict['pam'] = sequence[pam:pam+pam_len]
        cut_dict['guide'] = sequence[pam-20:pam+pam_len]
        cut_dict['guide_loc'] = (region[0], cut_dict['abs_break'] - 17, cut_dict['abs_break'] + 3 + pam_len)
        cut_dict['region'] = region

        if strand == '+':
            cut_dict['strand'] = '+'
        elif strand == '-':
            cut_dict['strand'] = '-'

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
    
    # TODO: check if positions on - strand are correct
    cuts_pos = fill_dict(seq, start, pams, pam_len, max_flanking_length, (chromosome, start, end, '+'))
    cuts_neg = fill_dict(rev_seq, start, rev_pams, pam_len, max_flanking_length, (chromosome, start, end, '-'))
    cut_sites = cuts_pos + cuts_neg
    
    return cut_sites


def parse_args():
    parser = argparse.ArgumentParser(description='Microhomology predictor.')
    
    parser.add_argument('--sequence-file', '-i', dest='sequence', help='File with the target sequence (TXT or FASTA).')
    parser.add_argument('--region', '-r', dest='region', help='Region in AgamP4 genome [2L:1530-1590].')
    parser.add_argument('--gene', '-G', dest='gene', help='Genome of interest (AgamP4.7 geneset).')
    parser.add_argument('--variants', '-v', dest='variants', help='VCF file with variants.')
    parser.add_argument('--pam', '-P', dest='pam', help='Protospacer adjacent motif (IUPAC format)', default='NGG')
    parser.add_argument('--threads', '-t', dest='n_threads', type=int, help='Number of threads used.', default=1)
    parser.add_argument('--max-flanking', '-M', type=int, dest='max_flanking_length', help='Max length of flanking region.', default=40)
    parser.add_argument('--min-flanking', '-m', type=int, dest='min_flanking_length', help='Min length of flanking region.', default=25)
    parser.add_argument('--length-weight', '-w', type=float, dest='length_weight', help='Length weight - used in scoring.', default=20.0)
    parser.add_argument('--max-offtargets', type=int, dest='max_offtargets', help='Max number of reported offtargets', default=100)
    parser.add_argument('--n-patterns', '-p', type=int, dest='n_patterns', help='Number of MH patterns used in guide evaluation.', default=5)
    parser.add_argument('--output-folder', '-o', dest='output_folder', help="Output folder.")
    parser.add_argument('--feature-type', '-f', dest='feature', help='Type of genomic feature to focus guide search on.', default=None)

    return parser.parse_args()


def define_genomic_region(chromosome, start, end):
    records = SeqIO.parse(os.path.join(ROOT_PATH, 'data', 'references', 'AgamP4.fa'), "fasta")
    reference = {}

    for record in records:
        reference[record.id] = record

    chr_seq = str(reference[chromosome].seq.upper())
    region = (chromosome, start, end, chr_seq)

    return region

def annotate_guides(cut_sites, ann_db, feature):
    '''
    Use GFF annotation to annotate guides
    (Optional - feature) Only keep guides that are located on a given type of genomic feature
    '''

    for cut_site in cut_sites[:]:
        location = cut_site['guide_loc']

        if feature == 'intergenic' or feature == 'intron':
            feature_types = [f[2] for f in ann_db.region(region=location)]

            if feature == 'intergenic' and 'gene' in feature_types:
                cut_sites.remove(cut_site)

            elif feature == 'intron' and 'gene' not in feature_types or feature == 'intron' and 'exon' in feature_types:
                cut_sites.remove(cut_site)

            else:
                features = [f for f in ann_db.region(region=location)]
                cut_site.update({'annotation': features})

        else:
            features = [f for f in ann_db.region(region=location, featuretype=feature)]
            if not features:
                cut_sites.remove(cut_site)
            cut_site.update({'annotation': features})

    return cut_sites


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------
def main():
    ascii_header = r'''
                                                                
                    ||||||            ||        ||            
                  ||        ||    ||        ||||||    ||||    
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
    feature = args.feature

    # ------------------------------------------------------------
    # Handle input arguments
    # ------------------------------------------------------------

    if args.region or args.gene:
        ann_db = gffutils.FeatureDB(os.path.join(ROOT_PATH, 'data', 'references', 'AgamP4.7'), keep_order=True)
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
            gene = ann_db[args.gene]
        except:
            logger.error('Gene not found: {}'.format(args.gene))
            quit()

        chromosome = gene.seqid
        start = gene.start
        end = gene.end

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

    # variants input -------------------------------------------------------
    if all(region) and args.variants:
        '''
        Option -v: get variants from VCF file
        '''
        with open(args.variants, 'r') as vcf_file:
            vcf_reader = vcf.Reader(vcf_file)

            logger.info('Reading VCF for region {}:{}-{}'.format(chromosome, start, end))
            variants = [v for v in vcf_reader.fetch(chromosome, start, end)]
    else:
        variants = []

    if not args.region and not args.sequence and not args.gene:
        logger.error('Please define the region of interest (-r) or provide the sequence (-i). Use -h for help.')
        quit()

    if args.feature is not None and ann_db.count_features_of_type(feature) == 0:
        feature_types = [f for f in ann_db.featuretypes()]
        logger.error('No feature of this type detected. Features present in current database are the following: {}.'.format(', '.join(feature_types)))
        quit()

    if ann_db:
        cut_sites = annotate_guides(cut_sites, ann_db, feature)

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

        logger.info('Simulating MMEJ ...')
        iterable_cut_sites = [(cut_site, args.n_patterns) for cut_site in cut_sites]
        cut_sites = list(tqdm(pool.istarmap(simulate_end_joining, iterable_cut_sites), total=len(iterable_cut_sites), ncols=100))
        
        logger.info('Add conservation and variation score ...')
        cut_sites = apply_conservation_variation_score(cut_sites, args.conservation_store, variants, pool)

        logger.info('Finding offtargets ...')
        targets_df = run_bowtie(cut_sites, args.max_offtargets, args.n_threads)

        pool.close()

        # create output dir if it doesn't exist
        # # create output dir if it doesn't exist
        # if not os.path.exists(args.output_folder):
        #     os.makedirs(args.output_folder)

        if args.sequence:
            # simple output
            save_guides_list_simple(cut_sites, args.output_folder, args.n_patterns)
            save_detailed_list_simple(cut_sites, args.output_folder, args.n_patterns)
        else:
            save_guides_list(cut_sites, args.output_folder, args.n_patterns)
            save_detailed_list(cut_sites, args.output_folder, args.n_patterns)
            save_to_bed(cut_sites, args.output_folder, args.n_patterns)

    else:
        logger.error('No output folder selected. Please define it by using -o option.')
        quit()
