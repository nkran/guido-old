import sys
import os
import requests
import json

import log

logger = log.createCustomLogger('requests')

SERVER = "https://www.vectorbase.org/rest"
HEADERS = {
    'text': 'text/plain',
    'json': 'application/json',
    'gff': 'text/x-gff3'
}


def parse_region(region):
    chromosome = region.split(':')[0]
    start = int(region.split(':')[1].split('-')[0])
    end = int(region.split(':')[1].split('-')[1])

    return (chromosome, start, end)


def make_get_request(path, header_type):

    path = os.path.join(SERVER, path)
    r = requests.get(path, headers={ "Content-Type" : HEADERS[header_type]})

    if not r.ok:
        r.raise_for_status()
        logger.error("Request failed.")
        sys.exit()

    if header_type == 'json':
        return json.loads(r.text)
    else:
        return r.text


def make_post_request(path, data, header_type):

    path = os.path.join(SERVER, path)
    r = requests.post(path, headers={"Content-Type" : HEADERS[header_type], "Accept": HEADERS[header_type]},
                            data=json.dumps(data))

    if not r.ok:
        r.raise_for_status()
        logger.error("Request failed.")
        sys.exit()

    if header_type == 'json':
        return json.loads(r.text)
    else:
        return r.text


def request_region_sequence(species, region):

    path = os.path.join('sequence', 'region', species, region)
    sequence = make_get_request(path, 'text')
    chromosome, start, end = parse_region(region)

    return (chromosome, start, end, sequence)


def request_gene_sequence(species, gene):

    path_seq = os.path.join('sequence', 'id', gene)
    path_info = os.path.join('overlap', 'id', gene + '?feature=gene')
    sequence = make_get_request(path_seq, 'text')
    gene_info = make_get_request(path_info, 'json')[0]

    return (gene_info['seq_region_name'], gene_info['start'], gene_info['end'], sequence)


def request_region(species, region, ftype):
    path = os.path.join('overlap', 'region', species, region + '?feature=' + ftype)
    info = make_get_request(path, 'json')

    return info  


def request_feature(species, gene, ftype):
    """ 
    Available types: Enum(gene, transcript, cds, exon, repeat, simple, misc, 
    variation, somatic_variation, structural_variation, 
    somatic_structural_variation, constrained, regulatory, segmentation, motif, 
    chipseq, array_probe)
    """
    path = os.path.join('overlap', 'id', '{}?feature={}'.format(gene, ftype))
    info = make_get_request(path, 'json')
    return info


def request_var_info(species, var_id):
    path = os.path.join('variation', species, var_id + '?pops=1')
    info = make_get_request(path, 'json')

    return info


def request_var_info_list(species, var_id_list):
    path = os.path.join('variation', species + '?pops=1')
    data = {'ids': var_id_list}
    info = make_post_request(path, data, 'json')
    
    return info