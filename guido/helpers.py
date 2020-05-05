import re
import argparse
import pandas

from collections import namedtuple

Region = namedtuple('Region', ['chromosome', 'start', 'end', 'sequence', 'annotation'])


def parse_args():
    parser = argparse.ArgumentParser(description='Microhomology predictor.')

    parser.add_argument(
        '--sequence-file',
        '-i',
        dest='sequence',
        help='File with the target sequence (TXT or FASTA).',
    )
    parser.add_argument(
        '--region',
        '-r',
        dest='region',
        help='Region in AgamP4 genome [2L:1530-1590].'
    )
    parser.add_argument(
        '--gene',
        '-G',
        dest='gene',
        help='Genome of interest (AgamP4.7 geneset).'
    )
    parser.add_argument(
        '--variants',
        '-v',
        dest='variation_store',
        help='VCF file with variants.',
        default=False,
    )
    parser.add_argument(
        '--conservation',
        '-c',
        dest='conservation_store',
        help='Path to Zarr store with conservation data.',
        default=False,
    )
    parser.add_argument(
        '--pam',
        '-P',
        dest='pam',
        help='Protospacer adjacent motif (IUPAC format)',
        default='NGG',
    )
    parser.add_argument(
        '--threads',
        '-t',
        dest='n_threads',
        type=int,
        help='Number of threads used.',
        default=1,
    )
    parser.add_argument(
        '--max-flanking',
        '-M',
        type=int,
        dest='max_flanking_length',
        help='Max length of flanking region.',
        default=40,
    )
    parser.add_argument(
        '--min-flanking',
        '-m',
        type=int,
        dest='min_flanking_length',
        help='Min length of flanking region.',
        default=25,
    )
    parser.add_argument(
        '--length-weight',
        '-w',
        type=float,
        dest='length_weight',
        help='Length weight - used in scoring.',
        default=20.0,
    )
    parser.add_argument(
        '--max-offtargets',
        type=int,
        dest='max_offtargets',
        help='Max number of reported offtargets',
        default=100,
    )
    parser.add_argument(
        '--n-patterns',
        '-p',
        type=int,
        dest='n_patterns',
        help='Number of MH patterns used in guide evaluation.',
        default=5,
    )
    parser.add_argument(
        '--disable-mmej',
        dest='disable_mmej',
        type=bool,
        nargs='?',
        const=True,
        help="Disable MMEJ prediction.",
        default=False,
    )
    parser.add_argument(
        '--disable-off-targets',
        dest='disable_offtargets',
        type=bool,
        nargs='?',
        const=True,
        help="Disable off-targets search.",
        default=False,
    )
    parser.add_argument(
        '--output-folder',
        '-o',
        dest='output_folder',
        help="Output folder."
    )
    parser.add_argument(
        '--feature-type',
        '-f',
        dest='feature',
        help='Type of genomic feature to focus guide search on.',
        default=None,
    )
    parser.add_argument(
        '--dump',
        dest='dump',
        help="Dump pickled cut_sites object to the output folder.",
        default=False,
    )
    parser.add_argument(
        '--genome-info',
        '-g',
        dest='genome_info',
        help='Genome information pickle file created using guido-build function.'
    )

    return parser.parse_args()


def chunks(l, n):
    """
    chunks a list into n-sized pieces
    """
    for i in range(0, len(l), n):
        # Create an index range for l of n items:
        yield l[i:i+n]


def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in seq[::-1]])


def parse_MD_tag(sequence, md_tag):
    tag = md_tag.split(':')[2]
    tag_chain = re.findall("([0-9]+)([AaCcGgTt]?)", tag)

    string = ''

    for link in tag_chain:
        dist, base = link

        string += '.' * int(dist)
        string += base

    return string


def parse_gff_info(feature, info):
    if feature in ['mRNA', 'CDS', 'exon', 'transcript']:
        props = info.replace(' ','').split(';')
        try:
            # GFF3 annotation attribute split
            info_dict = dict(s.split('=') for s in props)
        except:
            # GTF annotation attribute split
            info_dict = dict(s[:-1].split('"') for s in props if len(s) > 0)
    else:
        info_dict = {}

    return feature, info_dict
