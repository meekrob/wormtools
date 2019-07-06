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

def expand_match(m): 
    return m.start(), m.end(), m.string[m.start():m.end()]

def format_motif(arg):
    # copy an IUPAC-specified motif into its reverse complement
    # and merge the two into a legal regular expression
    usr_motif = arg.upper()
    
    # divide motif into a list of (symbol,range). A range is like '{0,3}'
    motif = []
    rangestr = ''

    for c in usr_motif:
        # add supported symbols
        #if c in 'ACGTMKWSYRN':
        if c in 'ACGTMKWSYRNBDHV':
            if rangestr:
                lc,r = motif[-1]
                motif[-1] = (lc,rangestr)
                rangestr = ''
            motif.append((c,None))
        # accumulate supported range character
        elif c in '{},0123456789':
            rangestr += c
        else:
            raise Exception("character " + c + " not a valid motif character") 

    # do the last one from the above loop 
    if rangestr:
        lc,r = motif_noranges[-1]
        motif_noranges[-1] = (lc,rangestr)
    
    # create the forward regular expression
    regex_motif = ''
    for char,ranges in motif:
        regex_motif += CODES[char]
        if ranges is not None:
            regex_motif += ranges
        
    #regex_motif += '|' # 'or'
    revcmp_motif = ''

    # create and add the reverse complemented regular expression
    motif.reverse()
    for char,ranges in motif:
        revcmp_motif += CODES[ char.translate(MOTIF_COMPLEMENT) ] 
        if ranges is not None:
            revcmp_motif += ranges

    return regex_motif,revcmp_motif

forward_pattern, reverse_pattern = format_motif( sys.argv[1] )
print >>sys.stderr, "searching", forward_pattern, "|", reverse_pattern

# prepare to read from file or stdin
mfastas = {}
current_fa = None
if len(sys.argv) > 2:
    inseq = open(sys.argv[2]) # must be single, not multi fasta

    print >>sys.stderr, "reading in sequences:",
    for line in inseq:
        if line.startswith('>'):
            current_fa = line.lstrip('>').strip()
        else:
            if current_fa not in mfastas: 
                mfastas[current_fa] = ''
            mfastas[current_fa] += line.strip().upper()

    print >>sys.stderr, "done."
        
else: 
    print >>sys.stderr, "%s motif file.fa" % sys.argv[0]
    sys.exit(0)

f_compiled = re.compile(forward_pattern)
r_compiled = re.compile(reverse_pattern)

print >>sys.stderr, "scanning sequences:",
for header,seq in mfastas.items():
    
    # search sequences
    f_result = [ expand_match(m) + ('+',) for m in f_compiled.finditer(seq)]
    r_result = [ expand_match(m) + ('-',) for m in r_compiled.finditer(seq)]
    result = f_result + r_result
    result.sort(key=lambda a: (a[0],a[1]))

    for i,r in enumerate(result): 
        s,e,seq,strand = r[0], r[1], r[2], r[3]
        print header, s,e, strand, seq

print >>sys.stderr, "done"
