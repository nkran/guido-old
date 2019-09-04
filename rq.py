from  guido.request import request_region_sequence, request_gene_sequence, request_variations, request_feature

# region = request_region_sequence('anopheles_gambiae', '2L:2422152-2422852')
# gene = request_gene_sequence('anopheles_gambiae', 'AGAP007280')

var = request_feature('anopheles_gambiae', 'AGAP004707?feature=variation')


# print type(var)

for v in var:
    if v['source'] == 'VBP0000227':
        print v
