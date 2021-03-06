import sys
"""
Static methods for parsing a line in GFF format.
============================================================================
Parse according to: http://www.ensembl.org/info/website/upload/gff.html

Fields must be tab-separated. 
Also, all but the final field in each feature line must contain a value; "empty" columns should be denoted with a '.'
----------------------------------------------------------------------------------------------------------------------
seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. 
source  - name of the program that generated this feature, or the data source (database or project name)
feature - feature type name, e.g. Gene, Variation, Similarity
start   - Start position of the feature, with sequence numbering starting at 1.
end     - End position of the feature, with sequence numbering starting at 1.
score   - A floating point value.
strand  - defined as + (forward) or - (reverse).
frame   - One of '0', '1' or '2'. 
attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.
----------------------------------------------------------------------------------------------------------------------

"""
# static data for tests
test_line_0 = 'IV	WormBase	gene	695	14926	.	+	.	gene_id "WBGene00021406"; gene_source "WormBase"; gene_biotype "protein_coding";\n'

# Return a dict keyed on attribute name. Values are always arrays.
def parseAttr(attr_str):
    """
    Parse an attribute field as read in GFF format. 
    This is relatively expensive, so be wary of calling by default.
    >>> parseLine(test_line_0)
    {'seqname': 'IV', 'end': 14926, 'name': 'gene', 'start': 695, 'frame': '.', 'source': 'WormBase', 'score': None, 'attr_str': 'gene_id "WBGene00021406"; gene_source "WormBase"; gene_biotype "protein_coding";', 'strand': '+'}
    """
    d = {}
    fields = [s.strip() for s in attr_str.split(';')]
    for i,f in enumerate(fields):
        if not f: continue # sometimes trailing semicolons can cause empty fields
        try:
            key,val = [s.strip('"') for s in f.strip().split()]
        except:
            print >>sys.stderr, "failed to parse field %d `%s' in fields [%s]" % (i,f,attr_str)
            raise
        if key not in d: d[key] = []
        d[key].append(val)
    
    return d

# Parse a qualified GFF line into a dict, optionally returning parsed dict of attributes.
def parseLine(gff_line, parse_attributes=False, exclude_features=[]):
    """
    Create a GFF dict from a line of GFF format.
    By default, {'attr_str': '...'} is returned instead of parsed {'attr':{...}} for performance reasons. 
    GFF.get_attr will call GFFParser.parseAttr if needed.
    >>> parseLine(test_line_0)
    {'seqname': 'IV', 'end': 14926, 'name': 'gene', 'start': 695, 'frame': '.', 'source': 'WormBase', 'score': None, 'attr_str': 'gene_id "WBGene00021406"; gene_source "WormBase"; gene_biotype "protein_coding";', 'strand': '+'}
    """
    gff = {}
    fields = gff_line.strip().split("\t")
    # sequence/source/feature
    try:
        gff['seqname'] = fields[0] # name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix.
        gff['source']  = fields[1] # name of the program the generated this record
        gff['name']    = fields[2] # feature type name, e.g. Gene, Variation, Similarity
    except IndexError:
        print >>sys.stderr, "not enough fields(%d) in line `%s'" % (len(fields),gff_line)
        raise
    if exclude_features and gff['name'] in exclude_features:
        return False
    # position in sequence
    gff['start']  = int(fields[3]) # Start position of the feature, with sequence numbering starting at 1.
    gff['end']    = int(fields[4]) # End position of the feature, with sequence numbering starting at 1.
    # score or Nonetype if '.'
    gff['score'] = None if fields[5] == '.' else float(fields[5])
    # strand/frame
    gff['strand'] = fields[6] # defined as + (forward) or - (reverse).
    gff['frame']  = fields[7] # One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
    # attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.
    if parse_attributes: # TODO: As this is the most expensive routine, make this false by default, and have the @staticmethod parse when necessary
        gff['attr'] = GFFParser.parseAttr(fields[8])
    else:
        gff['attr_str'] = fields[8]

    return gff

if __name__ == "__main__":
    import doctest
    doctest.testmod()
