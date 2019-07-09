from __future__ import print_function
import os
import log
import json

from collections import Counter
from request import request_var_info_list

logger = log.createCustomLogger('output')

# VARIATION_DICT = {'VBP0000227': 'G3',
#                   'VBP0000163': 'Ag1000g',
#                   'VBP0000002': 'AgSNP01_Bamako'}

def save_guides_list(cut_sites, output_folder, n_patterns):
    '''
    Print a list of guides
    '''

    filename = os.path.join(output_folder, 'guides_list_' + str(n_patterns) + '.txt')

    with open(filename, 'w') as f:
        logger.info('Creating table with guides: {}'.format(filename))

        print("guide_sequence\tgenomic_location\texon_name\tstrand\toff_target_analysis\tMMEJ_score\tMMEJ_sum_score\tMMEJ_top_score\tSNP_count\twt_prob\tSNP_info\tMMEJ_out_of_frame_del", file = f)
        
        for c in cut_sites:
            if 'complete_score' not in c.keys():
                print(c['seq'], 'N' in c['seq'])

        for guide_dict in sorted(cut_sites, key=lambda x: (x['complete_score'], x['sum_score']), reverse=True):
            mmej_frames = ('\t'.join('{}'.format(pattern['frame_shift']) for pattern in guide_dict['top_patterns']))
            location = guide_dict['guide_loc'][0] + ':' + str(guide_dict['guide_loc'][1]) + '-' + str(guide_dict['guide_loc'][2])

            variants_count = Counter([v['source'] for v in guide_dict['variants']])
            variants_string = " ".join(["{}:{}_".format(k, v) for k, v in variants_count.items()])[:-1]

            if 'annotation' in guide_dict.keys():
                exons = [e['Name'] for e in guide_dict['annotation'] if e.featuretype == 'exon']

                if exons:
                    exon_names = ",".join(exons[0])
                else:
                    exon_names = ""
            else:
                exon_names = ""

            if 'off_targets' not in guide_dict.keys():
                continue
            else:
                if guide_dict['top_patterns']:
                    print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(guide_dict['guide'],
                        location,
                        exon_names,
                        guide_dict['strand'],
                        guide_dict['off_targets']['count'],
                        guide_dict['complete_score'],
                        guide_dict['sum_score'],
                        guide_dict['top_patterns'][0]['pattern_score'],
                        variants_string,
                        mmej_frames), file = f)


def prepare_variants_block(variants):

    variants = [v['id'] for v in variants]
    # print(variants)

    a = request_var_info_list('anopheles_gambiae', variants)
    print(a)
        # variant_details = request_var_info('anopheles_gambiae', v['id'])

        # for pop in variant_details['populations']:
        #     print(pop)


def save_detailed_list(cut_sites, output_folder, n_patterns):
    '''
    Print detailed list of guides
    '''

    filename = os.path.join(output_folder, 'guides_list_details_' + str(n_patterns) + '.txt')

    with open(filename, 'w') as f:

        logger.info('Creating detailed table with MH patterns: {}'.format(filename))

        for guide_dict in sorted(cut_sites, key=lambda x: (x['complete_score'], x['sum_score']), reverse=True):
            location = guide_dict['guide_loc'][0] + ':' + str(guide_dict['guide_loc'][1]) + '-' + str(guide_dict['guide_loc'][2])

            print("Guide: {}\tLocation: {}\tStrand: {}".format(guide_dict['guide'], location, guide_dict['strand']), file = f)
            print("MMEJ score: {}\tMMEJ sum score: {}\tMMEJ top score: {}\n".format(guide_dict['complete_score'], guide_dict['sum_score'], guide_dict['top_patterns'][0]['pattern_score']), file = f)

            print("Top MMEJ patterns", file = f)
            print("Pattern\tScore\tDeletion size\tProduces out-of-frame deletion\tMH seq\tDeletion seq", file = f)
            for pattern in guide_dict['top_patterns']:
                print("{}\t{}\t{}\t{}\t{}\t{}".format(pattern['left'] + pattern['right'], pattern['pattern_score'], len(pattern['deletion_seq']), pattern['frame_shift'], pattern['pattern'], pattern['deletion_seq']), file = f)

            print("\nVariants", file = f)
            if guide_dict['variants'] and guide_dict['guide'] == 'GAGGAAGAAAGTGAGGAGGAGGG':
                guide_dict['annotation'] = {}
                print(guide_dict.keys())
                with open('/flash/repos/testing/b/7280_guide.json', 'w') as outfile:
                    json.dump(guide_dict, outfile)

                print("\nWT allele probability: {}".format(round(guide_dict['wt_prob'], 4)), file = f)
            else:
                print("No SNPs were found for given gRNA.", file = f)

            print("\nOff-targets", file = f)

            if sum(guide_dict['off_targets']['count']) > 0:

                print("Number off-targets with", file = f)
                print("[0, 1, 2, 3] mismatches", file = f)
                print(str(guide_dict['off_targets']['count']), file = f)

                print("\nChromosome\tStart position\tStrand\tMismatches", file = f)

                for ot in guide_dict['off_targets']['ot']:
                    mismatches = ot['mismatches'].replace(',', ', ').replace(':', ': ')
                    print("{}:{} {}\t{}".format(ot['chromosome'], ot['start'], ot['strand'], mismatches), file = f)
            else:
                print("No off-targets were found for given gRNA.", file = f)

            print("\n....................................................................................................", file = f)

# ---------------------
# simple output
# ---------------------

def save_guides_list_simple(cut_sites, output_folder, n_patterns):
    '''
    Print a simple list of guides
    '''

    filename = os.path.join(output_folder, 'guides_list_' + str(n_patterns) + '.txt')

    with open(filename, 'w') as f:

        logger.info('Creating table with guides: {}'.format(filename))

        print("guide_sequence\tgenomic_location\tstrand\toff_target_analysis\tMMEJ_score\tMMEJ_sum_score\tMMEJ_top_score\tMMEJ_out_of_frame_del", file = f)

        for guide_dict in sorted(cut_sites, key=lambda x: (x['complete_score'], x['sum_score']), reverse=True):
            mmej_frames = ('\t'.join('{}'.format(pattern['frame_shift']) for pattern in guide_dict['top_patterns']))
            location = guide_dict['guide_loc'][0] + ':' + str(guide_dict['guide_loc'][1]) + '-' + str(guide_dict['guide_loc'][2])

            if not 'off_targets' in guide_dict.keys():
                off_targets = 0
            else:
                off_targets = guide_dict['off_targets']['count']

            if guide_dict['top_patterns']:
                print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(guide_dict['guide'],
                    location,
                    guide_dict['strand'],
                    off_targets,
                    guide_dict['complete_score'],
                    guide_dict['sum_score'],
                    guide_dict['top_patterns'][0]['pattern_score'],
                    mmej_frames), file = f)


def save_detailed_list_simple(cut_sites, output_folder, n_patterns):
    '''
    Print detailed list of guides
    '''

    filename = os.path.join(output_folder, 'guides_list_details_' + str(n_patterns) + '.txt')

    with open(filename, 'w') as f:

        logger.info('Creating detailed table with MH patterns: {}'.format(filename))

        for guide_dict in sorted(cut_sites, key=lambda x: (x['complete_score'], x['sum_score']), reverse=True):
            location = guide_dict['guide_loc'][0] + ':' + str(guide_dict['guide_loc'][1]) + '-' + str(guide_dict['guide_loc'][2])

            print("Guide: {}\tLocation: {}\tStrand: {}".format(guide_dict['guide'], location, guide_dict['strand']), file = f)
            print("MMEJ score: {}\tMMEJ sum score: {}\tMMEJ top score: {}\n".format(guide_dict['complete_score'], guide_dict['sum_score'], guide_dict['top_patterns'][0]['pattern_score']), file = f)

            print("Top MMEJ patterns", file = f)
            print("Pattern\tScore\tDeletion size\tProduces out-of-frame deletion\tMH seq\tDeletion seq", file = f)
            for pattern in guide_dict['top_patterns']:
                print("{}\t{}\t{}\t{}\t{}\t{}".format(pattern['left'] + pattern['right'], pattern['pattern_score'], len(pattern['deletion_seq']), pattern['frame_shift'], pattern['pattern'], pattern['deletion_seq']), file = f)

            print("\nOff-targets", file = f)

            if 'off_targets' in guide_dict.keys() and sum(guide_dict['off_targets']['count']) > 0:

                print("Number off-targets with", file = f)
                print("[0, 1, 2, 3] mismatches", file = f)
                print(str(guide_dict['off_targets']['count']), file = f)

                print("\nChromosome\tStart position\tStrand\tMismatches", file = f)

                for ot in guide_dict['off_targets']['ot']:
                    mismatches = ot['mismatches'].replace(',', ', ').replace(':', ': ')
                    print("{}:{} {}\t{}".format(ot['chromosome'], ot['start'], ot['strand'], mismatches), file = f)
            else:
                print("No off-targets were found for given gRNA.", file = f)

            print("\n....................................................................................................", file = f)


def save_to_bed(cut_sites, output_folder, n_patterns):
    '''
    Print a list of guides
    '''

    filename = os.path.join(output_folder, 'guides_list_' + str(n_patterns) + '.bed')

    with open(filename, 'w') as f:

        logger.info('Creating table with guides: {}'.format(filename))

        print("guide_sequence\tgenomic_location\texon_name\tstrand\toff_target_analysis\tMMEJ_score\tMMEJ_sum_score\tMMEJ_top_score\tSNP_count\twt_prob\tSNP_info\tMMEJ_out_of_frame_del", file = f)

        for guide_dict in sorted(cut_sites, key=lambda x: (x['complete_score'], x['sum_score']), reverse=True):
            mmej_frames = ('\t'.join('{}'.format(pattern['frame_shift']) for pattern in guide_dict['top_patterns']))
            location = guide_dict['guide_loc'][0] + ':' + str(guide_dict['guide_loc'][1]) + '-' + str(guide_dict['guide_loc'][2])

            variants_count = len(guide_dict['variants'])
            variants_string = " ".join(["{}:{}/{}()".format(v['start'], v['alleles'][0], v['alleles'][0]) for v in guide_dict['variants']])

            if 'annotation' in guide_dict.keys():
                exons = [e['Name'] for e in guide_dict['annotation'] if e.featuretype == 'exon']

                if exons:
                    exon_names = ",".join(exons[0])
                else:
                    exon_names = ""
            else:
                exon_names = ""

            if guide_dict['top_patterns']:
                print("{}\t{}\t{}\t{}\t{}\t{}".format(
                    guide_dict['guide_loc'][0],
                    guide_dict['guide_loc'][1],
                    guide_dict['guide_loc'][2],
                    guide_dict['guide'],
                    guide_dict['complete_score'],
                    guide_dict['strand']), file = f)
