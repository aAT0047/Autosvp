import sys
import os
from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SubsMat import MatrixInfo

# the neoepitope file, one epitope sequence per line
file_epitope = sys.argv[1]

# the fasta file for epitopes in IEDB
file_iedb = sys.argv[2]

# outfile
outfile = sys.argv[3]

gap_open = -11
gap_extend = -1
matrix = MatrixInfo.blosum62

with open(file_epitope, 'r') as f:
    epitopes_neo = [Seq(line.rstrip(), IUPAC.protein) for line in f]

epitopes_iedb = list(SeqIO.parse(file_iedb, 'fasta'))

dic_neo = {}
for epitope_neo in epitopes_neo:
    maxscore = 0
    maxref = None
    maxlength = 0
    for epitope_iedb in epitopes_iedb:
        aln = pairwise2.align.localds(str(epitope_neo), str(epitope_iedb.seq), matrix, gap_open, gap_extend, one_alignment_only=True)
        if aln:
            aln = aln[0]
            if aln[2] > maxscore:
                maxscore = aln[2]
                maxref = epitope_iedb
                aln_detail = pairwise2.format_alignment(*aln).split('\n')
                maxlength = aln_detail[1].count('|')
            else:
                continue
    dic_neo[str(epitope_neo)] = [maxscore, maxref, maxlength]

with open(outfile, 'w') as f:
	f.write('\t'.join(['Neoepitope', 'Score', 'Length', 'IEDB_Sequence', 'Description']) + '\n')
	for epitope in dic_neo:
		if dic_neo[epitope][1]:
			f.write('\t'.join([epitope, str(dic_neo[epitope][0]), str(dic_neo[epitope][2]), str(dic_neo[epitope][1].seq), str(dic_neo[epitope][1].description)]) + '\n')
		else:
			f.write('\t'.join([epitope, str(dic_neo[epitope][0]), str(dic_neo[epitope][2]), 'None', 'None']) + '\n')
