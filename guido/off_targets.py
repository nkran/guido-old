import os
import subprocess

import log

logger = log.createCustomLogger('off-targets')
ROOT_PATH = os.path.abspath(os.path.dirname(__file__))

def run_bowtie(cut_sites, genome_index_path):

    logger.info('Running Bowtie for off-targets detection ...')

    missmatches = 3
    reported_aln = 10

    guides = ",".join([cut_site['guide'] for cut_site in cut_sites])

    # run bowtie alignment
    bowtie_command = 'bowtie -a -v {} -l 23 --suppress 6 {} -c {}'.format(missmatches, genome_index_path, guides)
    rproc = subprocess.Popen(bowtie_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)
    stdout, stderr = rproc.communicate()

    target_dict = {}

    for line in stdout.split(os.linesep):
        if line:
            b_id, b_strand, b_chromosome, b_start, b_seq, b_reps, b_mm = line.split('\t')

            if b_id not in target_dict.keys():
                target_dict[b_id] = []

            target = {
                'strand': b_strand,
                'chromosome': b_chromosome,
                'start': int(b_start),
                'reps': int(b_reps),
                'mismatches': b_mm
            }

            target_dict[b_id].append(target)

    with open(os.path.join(ROOT_PATH, 'data', 'bowtie_off_targets.txt'), 'w') as ot_file:
        ot_file.write(stdout)

    return target_dict


def off_target_evaluation(cut_sites, target_dict):

    logger.info('Evaluating off-targets ...')

    if not target_dict:
        for cut in cut_sites:
            cut.update({'off_targets': {
                        'count': [0, 0, 0, 0], # off-targets with 0, 1, 2, 3 mismatches
                        'ot': []
                        }})
    else:
        for tid in sorted(target_dict.keys(), key=lambda x: int(x)):

            targets = target_dict[tid]
            cut = cut_sites[int(tid)]

            cut_chromosome, cut_start, cut_end = cut['guide_loc']
            cut.update({'off_targets': {
                            'count': [0, 0, 0, 0], # off-targets with 0, 1, 2, 3 mismatches
                            'ot': []
                        }})

            for target in targets:
                # if the hit is not on-target
                if target['start'] not in range(cut_start - 1, cut_start + 1) and target['chromosome'] != cut_chromosome:

                    # append information about off-target
                    cut['off_targets']['ot'].append(target)

                    # split mismatch string by ',' and count the number of mismatches
                    # filter out empty string after splitting
                    missmatches = len(filter(None, target['mismatches'].split(',')))

                    # add to the mismatches counter
                    cut['off_targets']['count'][missmatches] += 1

    return cut_sites
