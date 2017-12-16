# convert GFF dict to array of BED6-compatible fields, including the [0,1]-based conversion
def GFF_as_BED6(gff):
    bed = [gff['seqname'], gff['start']-1,gff['end'],gff['name']]
    if gff['score'] is not None:
        bed.append( gff['score'] )
    else:
        bed.append( 0 )

    bed.append(gff['strand'])
    
    return bed

def gene_id(gff): # gff is normal dict
    return get_scalar_attr(gff, 'gene_id')
    
def get_scalar_attr(gff, attr_name): # gff is normal dict
    return get_attr(gff, attr_name)[0]

def get_attr(gff, attr_name):
    if not gff.has_key('attr'):
        if gff.has_key('attr_str'):
            gff['attr'] = GFFParser.parseAttr( gff['attr_str'] )
        else:
            raise Exception("gff has no 'attr' or 'attr_str' key. Was it produced via GFFParser?")

    return gff['attr'][attr_name]
