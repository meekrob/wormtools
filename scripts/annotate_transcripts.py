#!/usr/bin/env python
import GFF
import sys
from bx.bitset import BitSet # can be installed via pip
howmany = 0 # how many records to process, 0 for everything

#gene states
PROMOTER = "promoter"
INTRONIC = "intronic"
EXONIC = "exonic"
DOWNSTREAM_FLANK = "downstream_flank"

OUT_FH = open("genome.annotation.bed", "w")

#replacement for the following
# print "chr" + gene.seqname, s+base,e+base, gene.strand, current_state, gene.gene_id
def print_bed6(chrom, start, end, strand, label, gene_id):
    global OUT_FH
    score = 0
    name = gene_id + '.' + label
    print >>OUT_FH, "\t".join(map(str, [chrom, start, end, name, score, strand]))

def debug_print(*args):
    if howmany > 0:
        print >>sys.stderr, " ".join(map(str, args))

def debug_print_bitset(b, size):
    if howmany <= 0: return
    s = b.next_set(0)
    while s < size:
        e = b.next_clear(s)
        print >>sys.stderr, "debug_print_bitset", s, e
        if e == size: break
        s = b.next_set(e)

def createUnionBitset(a,b):
    union = a.clone()
    union.ior(b)
    return union

# subtract b from a, changing a in-place
def subtractBitset(a,b):
    b_comp = b.clone()
    b_comp.invert()
    a.iand(b_comp)

def to0baseEx(start1base, end1baseIncl):
    # when converting 1-based, inclusive such as GFF, 
    # to 0-base exclusive, such as BED
    return start1base-1, end1baseIncl

import pdb
from resolve_masks import three_states_two_masks # a tested submodule
def process_gene_model(gene, upstream_promoter_flank=2000, promoter_inset_flank=1000, downstream_flank=1000):
    if gene.isforward():
        if (gene.start - upstream_promoter_flank) < 1: # handle upstream_promoter_flank being less than left-end of chromosome
            upstream_promoter_flank = gene.start - 1
        gene_start,gene_end = to0baseEx(gene.start - upstream_promoter_flank,gene.end + downstream_flank)
    else:
        if (gene.start - downstream_flank) < 1: # handle downstream_flank being less than left-end of chromosome
            downstream_flank = gene.start - 1
        gene_start,gene_end = to0baseEx(gene.start - downstream_flank,gene.end + upstream_promoter_flank)

    base = gene_start
    size = gene_end - gene_start
    exon_mask = BitSet(size)
    promoter_mask = BitSet(size)
    #pdb.set_trace()
    # mark all exonic regions in all transcripts
    for t in gene.transcripts:
        exons = list(t.exons())
        for i_e, ex in enumerate(exons):
            # convert coordinates
            e_start,e_end = to0baseEx(ex['start'],ex['end'])
            pos = e_start - base
            l = e_end - e_start
            exon_mask.set_range(pos, l)
        if gene.isforward():
            promoter_exon = exons[0]
            p_start,p_end = to0baseEx(promoter_exon['start'],promoter_exon['end'])
            pos = p_start - base
            p_len = min(size, upstream_promoter_flank + promoter_inset_flank)
            promoter_mask.set_range( pos - upstream_promoter_flank, p_len )
        else:
            promoter_exon = exons[-1]
            p_start,p_end = to0baseEx(promoter_exon['start'],promoter_exon['end'])
            pos = p_end - base
            p_len = min(size, upstream_promoter_flank + promoter_inset_flank)
            try:
                promoter_mask.set_range( max(0, pos - promoter_inset_flank), p_len )
            except:
                print >>sys.stderr, "failed on gene", gene.gene_id, gene.start, gene.end, gene.strand, gene_start, gene_end, base, size
                raise

    things = three_states_two_masks(PROMOTER, promoter_mask, EXONIC, exon_mask, INTRONIC)

    # set the final 3' thing to be downstream flank
    threeprime_i = -1 if gene.isforward() else 0
    s,e,state = things[threeprime_i]
    if state == INTRONIC:
        things[threeprime_i] = (s,e,DOWNSTREAM_FLANK)

    for thing in things:
        s,e,current_state = thing
        print_bed6("chr" + gene.seqname, s+base,e+base, gene.strand, current_state, gene.gene_id)
        
def main():
    global OUT_FH
    infile = 'c_elegans.PRJNA13758.WS261.canonical_geneset.gtf'
    current_gene = None
    current_transcript = None
    transcript_count = 0
    for i,line in enumerate(open(infile)):
        # things to skip
        if line.startswith('#'): continue
        if line.startswith('MtDNA'): continue

        # refer to line if parsing fails
        try:
            gff = GFF.GFFParser.parseLine(line, exclude_features=['five_prime_utr', 'three_prime_utr', 'CDS','start_codon','stop_codon'])
        except:
            print >>sys.stderr, "failed on line %d:%s" % (i,line)
            raise
        if not gff: # False when feature is excluded (this allows the bypassing of parse_attr, which is expensive)
            continue

        # file is arranged by gene|transcript1[,transcript2,...]
        # so transcripts may be assumed to appear in groups
        if gff['name'] == 'transcript':
            #if current_transcript is not None: process_transcript(current_transcript)
            current_transcript = GFF.Transcript(gff)
            transcript_count += 1

            if current_gene is not None: 
                current_gene.add_transcript( current_transcript )

            # break if debugging
            if howmany and transcript_count > howmany: break
            
        elif gff['name'] == 'gene':
            if current_gene is not None:
                #print current_gene, len(current_gene)
                process_gene_model( current_gene )
                """
                for ptuple in list(current_gene.get_promoters(300)):
                    plist = list(ptuple)
                    plist[0] = "chr" + plist[0]
                    print " ".join(map(str,plist))
                """
            current_gene = GFF.Gene(gff)

        else: # add any other gff type (exon,CDS,etc) to current transcript
            current_transcript.append(gff)

    OUT_FH.close()

main()

def process_transcript(transcript):
    print transcript

