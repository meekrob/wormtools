#!/usr/bin/env python
import re,sys,string

if len(sys.argv) < 2:
    print sys.argv[0], "motif", "infile.fasta"
    print "    Convert motif into regular expression, plus reverse complement, and scan input sequence."
    print "If infile.fasta is not supplied, sequence is assumed to be from stdin"
    print "    motif can contain all IUPAC characters https://en.wikipedia.org/wiki/Nucleic_acid_notation"
    print "    A valid regular expression range, i.e. {min,max} may be specified to apply to the preceding character."
    print "    See https://www.gnu.org/software/grep/manual/grep.html#Fundamental-Structure for exact syntax."
    sys.exit(1)

VALID_CHARS = 'ACGTKRWMYSBDHVN'

CODES = { # use this to make a legal regex
    'A':'A', 'C':'C', 'G':'G', 'T':'T',
    'K': '[GT]', 'R': '[AG]', 'W': '[AT]',
    'M': '[AC]', 'Y': '[CT]', 'S': '[GC]', 
    'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]',
    'N': '[ACGT]' }

# translation tables for complementation
DNA_COMPLEMENT = string.maketrans('ACGT','TGCA')
MOTIF_COMPLEMENT = string.maketrans('ACGT' + 'KRWMYS' + 'BDHVN', 'TGCA' + 'MYWKRS' + 'VHDBN')

def format_motif(arg):
    usr_motif = arg.upper()
    
    # disassemble string into array with tuples (symbol, range) where range is a string like '{0,3}' that can be attached to a regex pattern
    motif = []
    rangestr = ''

    for c in usr_motif:
        if c in 'ACGTMKWSYRN':
            if rangestr:
                lc,r = motif[-1]
                motif[-1] = (lc,rangestr)
                rangestr = ''
            motif.append((c,None))
        elif c in '{},0123456789':
            rangestr += c
        else:
            raise Exception("character " + c + " not a valid motif character") 

    if rangestr:
        lc,r = motif_noranges[-1]
        motif_noranges[-1] = (lc,rangestr)

    
    regex_motif = ''
    for char,ranges in motif:
        regex_motif += CODES[char]
        if ranges is not None:
            regex_motif += ranges
        
    regex_motif += '|'
    motif.reverse()
    for char,ranges in motif:
        regex_motif += CODES[ char.translate(MOTIF_COMPLEMENT) ] 
        if ranges is not None:
            regex_motif += ranges

    return regex_motif

pattern = format_motif( sys.argv[1] )
print >>sys.stderr, "searching", pattern

# prepare to read from file or stdin
if len(sys.argv) > 2:
    inseq = open(sys.argv[2]) # must be single, not multi fasta
    header = fastafile.readline()
else: 
    inseq = sys.stdin

print >>sys.stderr, "reading in sequence:",
seq = ''
for line in inseq: seq += line.strip().upper()
# verify sequence
print >>sys.stderr, len(seq), "characters: `%s'" % seq
# search sequence
result = re.compile(pattern).finditer(seq)
for i,r in enumerate(result): 
    s,e = r.start(), r.end()
    print s,e, seq[s:e]

