from __future__ import print_function
import os
import log

logger = log.createCustomLogger('output')

def save_guides_list(cut_sites, output_folder, n_patterns):
    '''
    Print a list of guides
    '''

    filename = os.path.join(output_folder, 'guides_list_' + str(n_patterns) + '.txt')

    with open(filename, 'w') as f:

        logger.info('Creating table with guides: {}'.format(filename))

        print("guide_sequence\tgenomic_location\texon_name\tstrand\toff_target_analysis\tMMEJ_score\tMMEJ_sum_score\tMMEJ_top_score\tSNP_count\twt_prob\tMMEJ_out_of_frame_del", file = f)

        for guide_dict in sorted(cut_sites, key=lambda x: (x['complete_score'], x['sum_score']), reverse=True):
            mmej_frames = ('\t'.join('{}'.format(pattern['frame_shift']) for pattern in guide_dict['top_patterns']))
            location = guide_dict['guide_loc'][0] + ':' + str(guide_dict['guide_loc'][1]) + '-' + str(guide_dict['guide_loc'][2])

            variants_count = len(guide_dict['variants'])
            variants_string = " ".join(["{}:{}/{}({})".format(v.POS, v.REF, v.ALT, [round(v, 4) for v in v.aaf]) for v in guide_dict['variants']])

            exons = [e['Name'] for e in guide_dict['annotation'] if e.featuretype == 'exon']

            if exons:
                exon_names = ",".join(exons[0])
            else:
                exon_names = ""

            if guide_dict['top_patterns']:
                print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(guide_dict['guide'],
                    location,
                    exon_names,
                    guide_dict['strand'],
                    guide_dict['status'],
                    guide_dict['complete_score'],
                    guide_dict['sum_score'],
                    guide_dict['top_patterns'][0]['pattern_score'],
                    variants_count,
                    round(guide_dict['wt_prob'], 4),
                    variants_string,
                    mmej_frames), file = f)


# print the list of guides and their MH patterns
def save_mh_list(cut_sites, output_folder, n_patterns):
    '''
    Print the list of guides
    '''

    filename = os.path.join(output_folder, 'guides_list_mh_' + str(n_patterns) + '.txt')

    with open(filename, 'w') as f:

        logger.info('Creating detailed table with MH patterns: {}'.format(filename))

        for guide_dict in sorted(cut_sites, key=lambda x: (x['complete_score'], x['sum_score']), reverse=True):
            location = guide_dict['guide_loc'][0] + ':' + str(guide_dict['guide_loc'][1]) + '-' + str(guide_dict['guide_loc'][2])

            print("{}\t{}\t{}\t{}\t{}\t{}".format(guide_dict['guide'], location, guide_dict['strand'], guide_dict['status'], guide_dict['complete_score'], guide_dict['sum_score']), file = f)
            for pattern in guide_dict['top_patterns']:
                print("{}\t{}\t{}\t{}\t{}\t{}".format(pattern['left'] + pattern['right'], pattern['pattern_score'], len(pattern['deletion_seq']), pattern['frame_shift'], pattern['pattern'], pattern['deletion_seq']), file = f)

            print("....................................................................", file = f)
