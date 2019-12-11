import os
import argparse
import pickle
import subprocess

import guido.log as log

logger = log.createCustomLogger('root')


def parse_args():
    parser = argparse.ArgumentParser(description='Build genome index for bowtie to run alignment of off-targets.')

    parser.add_argument('--genome-file', '-g', dest='genome_file', help='Filename with the target sequence (.fa).')
    parser.add_argument('--genome-name', '-n', dest='genome_name', help='Genome base name')
    parser.add_argument('--annotation-file', '-a', dest='annotation_file', help='Filename with genome annotation (.gtf or .gff3).')

    return parser.parse_args()


def main():
    """
    TODO
    - handle only specific formats (.fa, .gz / .gtf, .gff3, .gz)
    - Convert fasta to twoBit compressed file
    - enable disabling (eg. twoBit)
    """

    ascii_header = r'''
                   ________  __________  ____        ____  __  ________    ____ 
                  / ____/ / / /  _/ __ \/ __ \      / __ )/ / / /  _/ /   / __ \
                 / / __/ / / // // / / / / / /_____/ __  / / / // // /   / / / /
                / /_/ / /_/ // // /_/ / /_/ /_____/ /_/ / /_/ // // /___/ /_/ / 
                \____/\____/___/_____/\____/     /_____/\____/___/_____/_____/  
                                                                                
                    '''

    print(ascii_header)

    logger.info("Building genome index so you don't have to.")
    args = parse_args()

    # Define genome_info dict structure for pickling
    genome_info = {'genome_name': '',
                   'description': input('Genome description:'),
                   'genome_file': '',
                   'annotation_file': ''}

    # ------------------------------------------------------------
    # Handle input arguments
    # ------------------------------------------------------------

    if not args.genome_file or not os.path.exists(args.genome_file):
        logger.error('Please use the -g argument to direct guido-build to the genome fasta file.')
        quit()
    else:
        genome_file_abspath = os.path.abspath(args.genome_file)
        genome_file_dirname = os.path.dirname(genome_file_abspath)
        genome_info['genome_file'] = genome_file_abspath

    if not args.genome_name:
        logger.error('Please use the -n argument to give your genome index files a name.')
        quit()
    else:
        genome_info['genome_name'] = args.genome_name

    if args.annotation_file:
        if not os.path.exists(args.annotation_file):
            logger.error('Annotation file does not exist. Check path is entered correctly.')
            quit()
        else:
            genome_info['annotation_file'] = os.path.abspath(args.annotation_file)

    # ------------------------------------------------------------
    # Execute
    # ------------------------------------------------------------

    # Create bowtie index
    bowtie_index_command = 'bowtie-build {} {}'.format(genome_file_abspath, os.path.join(genome_file_dirname, args.genome_name))
    rproc = subprocess.Popen(bowtie_index_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    stdout, stderr = rproc.communicate()

    # Create fai index
    # '''
    # Unnecessary for current usage.
    # May be useful when indexing larger projects.
    # '''
    # fai_index_command = 'samtools faidx {}'.format(args.genome_file)
    # print(fai_index_command.split())
    # rproc = subprocess.Popen(fai_index_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    # stdout, stderr = rproc.communicate()

    # Pickle genome_info dictionary
    with open(os.path.join(genome_file_dirname, 'genome_info.pickle'), 'wb') as f:
        pickle.dump(genome_info, f, protocol=pickle.HIGHEST_PROTOCOL)
