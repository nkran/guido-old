import argparse
import os
import pickle
import subprocess

from pyfaidx import Faidx

from guido.log import createCustomLogger

logger = createCustomLogger('root')


def parse_args():
    parser = argparse.ArgumentParser(description='Build genome index for bowtie to run alignment of off-targets.')

    parser.add_argument('--genome-file', '-g', dest='genome_file', help='Filename with the target sequence (FASTA).')
    parser.add_argument('--genome-name', '-n', dest='genome_name', help='Genome base name')
    parser.add_argument('--annotation-file', '-a', dest='annotation_file', help='Filename with genome annotation (.gtf, .gff3, .gz).')
    parser.add_argument('--threads', '-t', dest='n_threads', type=int, help='Number of threads used.', default=1)
    parser.add_argument('--description', '-d', dest='description', type=str, help='Description of genome.')

    return parser.parse_args()


def main():
    ascii_header = r'''
                   ________  __________  ____        ____  __  ________    ____ 
                  / ____/ / / /  _/ __ \/ __ \      / __ )/ / / /  _/ /   / __ \
                 / / __/ / / // // / / / / / /_____/ __  / / / // // /   / / / /
                / /_/ / /_/ // // /_/ / /_/ /_____/ /_/ / /_/ // // /___/ /_/ / 
                \____/\____/___/_____/\____/     /_____/\____/___/_____/_____/  
                                                                                
                    '''

    print(ascii_header)

    logger.info("Building genome indices so you don't have to. This will take a few minutes ... ")
    args = parse_args()

    # Define genome_info dict structure for pickling
    genome_info = {'genome_name': '',
                   'description': '',
                   'genome_file': '',
                   'annotation_file': '',
                   'fai_file': '',
                   'tbi_file': ''}

    # ------------------------------------------------------------
    # Handle input arguments
    # ------------------------------------------------------------

    if not args.genome_file:
        logger.error('Please use the -g argument to direct guido-build to the genome fasta file.')
        quit()
    else:
        if os.path.exists(args.genome_file):
            supported_seq_ext = ['.fa', '.fasta'] #TODO - check if the following formats are also accepted ['.fna', '.ffn', '.faa', '.frn']
            gen_seq_ext = os.path.splitext(args.genome_file)[1]
            if gen_seq_ext in supported_seq_ext:
                genome_file_abspath = os.path.abspath(args.genome_file)
                genome_file_dirname = os.path.dirname(genome_file_abspath)
                genome_info['genome_file'] = genome_file_abspath
            else:
                logger.error(f'Genome file format must be in a supported FASTA format, e.g. {" ".join(supported_seq_ext)}. Files must be unzipped in order to build bowtie indices.')
                quit()

    if not args.genome_name:
        logger.error('Please use the -n argument to give your genome index files a name.')
        quit()
    else:
        genome_info['genome_name'] = args.genome_name

    if args.annotation_file:
        ann_abspath = os.path.abspath(args.annotation_file)
        if os.path.exists(args.annotation_file):
            supported_ann_ext = ['.gtf', '.gff3'] # TODO - implement '.gz' format handling
            ann_path, ann_ext = os.path.splitext(ann_abspath)
            if ann_ext in supported_ann_ext:
                genome_info['annotation_file'] = ann_abspath
            else:
                logger.error('Annotation file format must be either in either GTF or GFF3 format.')
        else:
            logger.error('Annotation file does not exist. Check path is entered correctly.')
            quit()

    if args.description:
        genome_info['description'] = args.description

    if args.n_threads > 1:
        logger.info(f'Guido-build is dancing with {args.n_threads} threads ...')
    else:
        logger.info('Guido-build is dancing with a lonely single thread ...')

    # ------------------------------------------------------------
    # Execute
    # ------------------------------------------------------------

    # Create tabix files
    logger.info('Creating tabix files.')

    sorted_ann_abspath = ann_path + ".sorted" + ann_ext
    ann_sort_bgz_cmd = f'sort -k1,1 -k4,4n {ann_abspath} | bgzip > {sorted_ann_abspath}.gz'
    sort_proc = subprocess.run(ann_sort_bgz_cmd, shell=True)

    tabix_cmd = f'tabix {sorted_ann_abspath}.gz'
    tabix_proc = subprocess.run(tabix_cmd, stderr=subprocess.PIPE, shell=True)

    genome_info['tbi_file'] = f'{sorted_ann_abspath}.gz.tbi'

    # Create bowtie index
    logger.info('Building bowtie index.')
    bowtie_index_command = f'bowtie-build {genome_file_abspath} {os.path.join(genome_file_dirname, args.genome_name)} --threads {args.n_threads}'
    rproc = subprocess.Popen(bowtie_index_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    stdout, stderr = rproc.communicate()

    # Create fai index
    logger.info('Building fasta index.')
    Faidx(genome_file_abspath)
    genome_info['fai_file'] = f'{genome_file_abspath}.fai'

    # Pickle genome_info dictionary
    with open(os.path.join(genome_file_dirname, 'genome_info.pickle'), 'wb') as f:
        pickle.dump(genome_info, f, protocol=pickle.HIGHEST_PROTOCOL)
