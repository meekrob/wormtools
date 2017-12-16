@staticmethod
def GFF_as_BED6(gff):
    """
    Convert GFF dict to array of BED6-compatible fields, including the [0,1]-based conversion
    """
    bed = [gff['seqname'], gff['start']-1,gff['end'],gff['name']]
    if gff['score'] is not None:
        bed.append( gff['score'] )
    else:
        bed.append( 0 )

    bed.append(gff['strand'])
    
    return bed

def gene_id(gff): # gff is normal dict
    """
    Convenience method to get the scalar gene_id out of a GFF returned by GFF.GFFParser
    >>> gff = {'attr':{'gene_id':['sample_gene_id']}}
    >>> gene_id(gff)
    'sample_gene_id'
    """
    return get_scalar_attr(gff, 'gene_id')
    
def get_scalar_attr(gff, attr_name): # gff is normal dict
    """
    Convenience method to get a scalar attribute out of a GFF returned by GFF.GFFParser
    >>> gff = {'attr':{'sample_scalar_attr_name':['sample_scalar_attr_value']}}
    >>> get_scalar_attr(gff, 'sample_scalar_attr_name')
    'sample_scalar_attr_value'
    """
    return get_attr(gff, attr_name)[0]

def get_attr(gff, attr_name):
    """
    Convenience method to return an attribute list from a GFF returned by GFF.GFFParser, or parse the 'attr_str' if it hasn't happened yet.
    >>> gff = {'attr':{'sample_attr_name':['sample_attr_value']}}
    >>> get_attr(gff, 'sample_attr_name')
    ['sample_attr_value']
    """
    if not gff.has_key('attr'):
        if gff.has_key('attr_str'):
            gff['attr'] = GFFParser.parseAttr( gff['attr_str'] )
        else:
            raise Exception("gff has no 'attr' or 'attr_str' key. Was it produced via GFFParser?")

    return gff['attr'][attr_name]

if __name__ == "__main__":
    import doctest
    doctest.testmod()
