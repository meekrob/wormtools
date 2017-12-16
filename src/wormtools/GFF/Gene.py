"""
Accumulates GFF.Transcript objects and a GFF dict whose feature is 'gene', supporting direct access for seqname, start, end, gene_id, strand, and others.
"""
class Gene:
    def __init__(self,genegff):
        """
        Initialize on a GFF dict from GFF.GFFParseline whose name == 'gene'
        """
        if genegff['name'] != 'gene':
            raise Exception("constructor requires gff with name == 'gene'.")

        # things I care about at the moment
        self._seqname = genegff['seqname']
        self._start = genegff['start']
        self._end = genegff['end']
        self._gene_id = gene_id(genegff)
        self._strand =  genegff['strand']
        self._transcripts = [] 

    @property
    def transcripts(self):
        """List of transcripts vetted by add_transcript() """
        return self._transcripts
    @property
    def strand(self):
        """From genegff['strand'] """
        return self._strand
    @property
    def gene_id(self):
        """From genegff['gene_id'] """
        return self._gene_id
    @property
    def end(self):
        """From genegff['end'] """
        return self._end
    @property
    def start(self):
        """From genegff['start'] """
        return self._start
    @property
    def seqname(self): 
        """From genegff['seqname'] """
        return self._seqname

    def isforward(self):
        """True/False is strand '+'?"""
        return self.strand == '+'

    def add_transcript(self, transcript):
        """
        Add a GFF.Transcript object. 
        An unmatched gene_id between self and transcript throws AssertionError """
        assert self.gene_id == transcript.gene_id
        self.transcripts.append(transcript)

    def get_promoters(self, bp_offset):
        """
        Determine the 5' exon for each transcript and extend its start or end by bp_offset, 
        for '+' and '-' strand transcripts respectively.
        """
        for transcript in self.transcripts:
            exon = transcript.get_5prime_exon()
            if self.strand == '+':
                yield exon['seqname'], exon['start'] - bp_offset, exon['start'], "+promoter", self.gene_id + "|" + transcript.transcript_id
            else:
                yield exon['seqname'], exon['end'], exon['end'] + bp_offset, "-promoter", self.gene_id + "|" + transcript.transcript_id
            
    def get_transcript_ids(self):
        """See Transcript.transcript_id"""
        for transcript in self.transcripts:
            yield transcript.transcript_id

    def __len__(self):
        """Number of transcripts"""
        return len(self.transcripts)

    def __str__(self):
        return "<<<<<<<<<< GENE: %s >>>>>>>>>>>>" % self.gene_id

