import sys

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
    

class Gene:
    
    def __init__(self,genegff):
        if genegff['name'] != 'gene':
            raise Exception("constructor requires gff with name == 'gene'.")

        # things I care about at the moment
        self.seqname = genegff['seqname']
        self.start = genegff['start']
        self.end = genegff['end']
        self.gene_id = gene_id(genegff)
        self.strand =  genegff['strand']
        self.transcripts = [] 

    def isforward(self):
        return self.strand == '+'

    def add_transcript(self, transcript):
        assert self.gene_id == transcript.gene_id
        self.transcripts.append(transcript)

    def get_promoters(self, bp_threshold):
        for transcript in self.transcripts:
            exon = transcript.get_5prime_exon()
            if self.strand == '+':
                yield exon['seqname'], exon['start'] - bp_threshold, exon['start'], "+promoter", self.gene_id + "|" + transcript.transcript_id
            else:
                yield exon['seqname'], exon['end'], exon['end'] + bp_threshold, "-promoter", self.gene_id + "|" + transcript.transcript_id
            

    def get_transcript_ids(self):
        for transcript in self.transcripts:
            yield transcript.transcript_id

    def __len__(self):
        return len(self.transcripts)

    def __str__(self):
        return "<<<<<<<<<< GENE: %s >>>>>>>>>>>>" % self.gene_id

class Transcript:
    ALLOWED_FEATURES = ['CDS','exon','start_codon', 'stop_codon','five_prime_utr','three_prime_utr','transcript']

    def __init__(self,gffdata):
        self.gffdata = [] 
        self.gene_id = None
        self.transcript_id = None
        self.strand = None
        if type(gffdata) == type({}): # should be a single gff dict, as returned by a single call to GFFParser.parseLine
            self.append(gffdata)
        elif type(gffdata) == type([]): # should be a list of gff dicts
            self.extend(gffdata)
        else:
            raise Exception("Transcript was initialized by an unexpected type of data %s" % str(gffdata))

    def exons(self):
        for gff in self.gffdata:
            if gff['name'] == 'exon':
                yield gff

    def get_5prime_exon(self):
        return self.__get_exon_end(True)
    def get_3prime_exon(self):
        return self.__get_exon_end(False)

    def __get_exon_end(self, is_5prime):
        if self.strand is None:
            raise Exception("call to get_5prime_end() was made before being given a transcript gff or strand")

        exons = [e for e in self.exons()]
        
        if is_5prime:
            exons.sort(lambda a,b: cmp(a['start'],b['start']))
        else:
            exons.sort(lambda a,b: cmp(a['end'],b['end']))

        if self.strand == "+" and is_5prime:
            exons.sort(lambda a,b: cmp(a['start'],b['start']))
            return exons[0]
        elif self.strand == "+" and is_3prime:
            exons.sort(lambda a,b: cmp(a['end'],b['end']))
            return exons[-1]
        elif self.strand == "-" and is_5prime:
            exons.sort(lambda a,b: cmp(a['end'],b['end']))
            return exons[-1]
        else:
            exons.sort(lambda a,b: cmp(a['start'],b['start']))
            return exons[0]

    def append(self, gffdatum):
        if gffdatum['name'] == "transcript":
            self.transcript = gffdatum
            #self.gene_id = get_attr(self.transcript['attr']['gene_id'][0]
            self.gene_id = gene_id(self.transcript)

            #self.transcript_id = self.transcript['attr']['transcript_id'][0]
            self.transcript_id = get_scalar_attr(self.transcript, 'transcript_id') 
            self.strand = self.transcript['strand']
            self.start = self.transcript['start']
            self.end = self.transcript['end']
        elif gffdatum['name'] in self.ALLOWED_FEATURES:
            self.gffdata.append(gffdatum)
        else:
            raise Exception("feature type %s not allowed in transcript" % gffdatum['name'])

        self.gffdata.sort(lambda a,b: cmp(a['start'], b['start']))

    def extend(self, data):
        for datum in data: self.append(datum)

    def __str__(self):
        s = '%s:%s\n' % (self.transcript.gene_id, self.transcript.transcript_id)
        for gff in self.gffdata:
            s += "%s:%s %d %d %s\n" % (gff['seqname'], gff['name'].ljust(11), gff['start'], gff['end'], gff['strand'])

        return s

class GFFParser:
    """Parse according to: http://www.ensembl.org/info/website/upload/gff.html

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

    # Return a dict keyed on attribute name. Values are always arrays.
    @staticmethod
    def parseAttr(attr_str):
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
    @staticmethod
    def parseLine(gff_line, parse_attributes=False, exclude_features=[]):
        gff = {}
        fields = gff_line.strip().split("\t")
        # sequence/source/feature
        gff['seqname'] = fields[0] # name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix.
        gff['source']  = fields[1] # name of the program the generated this record
        gff['name']    = fields[2] # feature type name, e.g. Gene, Variation, Similarity
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

