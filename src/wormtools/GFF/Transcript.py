import sys
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



