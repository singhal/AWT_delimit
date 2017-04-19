import glob
import re
import os
import copy
import random

dir = '/Users/sonal/Desktop/dryad_34274/data/assembled/'
outdir = '/Users/Sonal/Desktop/Pseudemoia_entrecasteauxii/'

ambig = {	'M': ['A', 'C'], 'R': ['A', 'G'],
			'W': ['A', 'T'], 'S': ['C', 'G'],
			'Y': ['C', 'T'], 'K': ['G', 'T']
}

def get_seq(file):
	id = ''
	seq = {}
	o = open(file, 'r')
	for l in o:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			seq[id] = ''
		else:
			seq[id] += l.rstrip()
	o.close()
	return seq

pe = 'SP03_indexing12'
pen = 'PE70'

files = glob.glob(dir + '*fa')
for file in files:
	locname = re.search('(ENS.*)_ambig', file).group(1)
	seq = get_seq(file)
	if pe in seq:
		out = os.path.join(outdir, '%s.aln.fasta' % locname)
		o = open(out, 'w')
		seq1 = list(copy.copy(seq[pe]))
		seq2 = list(copy.copy(seq[pe]))

		for ix, base in enumerate(seq[pe]):
			if base not in ['A', 'T', 'C', 'G', 'N']:
				bases = copy.copy(ambig[base])
				random.shuffle(bases)

				seq1[ix] = bases[0]
				seq2[ix] = bases[1]

		o.write('>%s_1\n%s\n' % (pen, ''.join(seq1)))
		o.write('>%s_2\n%s\n' % (pen, ''.join(seq2)))

		o.close()