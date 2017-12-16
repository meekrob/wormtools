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

