import re
import glob
import os
import random
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outdir", required=True, help='Outdir.')
parser.add_argument("-r", "--run", required=True, help='Run number.')
parser.add_argument("-f", "--file", required=True, help="Fle to run.")
parser.add_argument("-m", "--miss", required=True, help="Missingness.")
args = parser.parse_args()

num_loci = 200
minpic = 0.001

def get_sp_loci(args):
        d = pd.read_csv(args.file)
        d = d.groupby('lineage')
        groups = {}
        for lineage, group in d:
                ind = random.sample(group['sample'].tolist(), 1)[0]
		ix = random.choice([1,2])
		groups['%s_%s' % (ind, ix)] = lineage
	
        outdir = args.outdir
	outdir = os.path.join(outdir, 'phylonet')
	if not os.path.isdir(outdir):
		os.mkdir(outdir)
	outdir = os.path.join(outdir, 'run%s' % args.run)
	if not os.path.isdir(outdir):
		os.mkdir(outdir)
	
        loc_file = os.path.join(args.outdir, 'locus_data.csv')
        d = pd.read_csv(loc_file)

        loci = d.ix[d.missingness >= float(args.miss), 'locus'].tolist()
       
        return outdir, groups, loci

def get_seq(aln, groups):
	f = open(aln, 'r')
	s = {}
	length = None

	for l in f:
		if not re.search('\s+\d+', l):
			(id, seq) = re.split('\s+', l.rstrip())
			s[id] = list(seq.upper())
			length = len(seq)
	f.close()

	keep = {}
	for id, seq in s.items():
		if id in groups:
			keep[id] = seq

	return keep, length					

def pic(bases):
	bases = [x for x in bases if x not in ['-', 'n', 'N']]
	
	counts = {}
	for i in list(set(bases)):
		count = bases.count(i)
		if count > 1:
			counts[i] = count

	if len(counts) > 1:
		return 1
	else:
		return 0


def get_pics(args, loci, groups):
	loc_data = {}
	for locus in loci:
		infile = os.path.join(args.outdir, 'alignments', '%s.aln.phy' % locus)
		seq, seqlen = get_seq(infile, groups)
		n_pics = 0
		for i in range(0, seqlen):
			bases = [seq[id][i] for id in seq]
			n_pics += pic(bases)
		loc_data[locus] = {'seq': seq, 'seqlen': seqlen, 
				   'pic': n_pics / float(seqlen), 'nind': len(seq)}

	return loc_data

outdir, groups, loci = get_sp_loci(args)
loc_data = get_pics(args, loci, groups)
# for locus in loc_data:
#	print(locus, loc_data[locus]['pic'], loc_data[locus]['seqlen'])
possible = []
for locus in loc_data:
	if loc_data[locus]['pic'] >= minpic and loc_data[locus]['seqlen'] > 150:
		possible.append(locus)
sel = random.sample(possible, num_loci)

tot_len = sum([loc_data[locus]['seqlen'] for loc in sel])
num_inds = len(groups)

out = os.path.join(outdir, 'phylonet_miss%s_pic%s.nex' % (args.miss, minpic))
o = open(out, 'w')
o.write('#NEXUS\n')
o.write('Begin data;\n')
o.write('\tDimensions ntax=%s nchar=%s;\n' % (num_inds, tot_len))
o.write('\tFormat datatype=dna symbols="ACTG" missing=? gap=-;\n')
o.write('\tMatrix\n')

for loc_name in sel:
	seq = loc_data[loc_name]['seq']
	length = loc_data[loc_name]['seqlen']
	o.write('[%s, %s]\n' % (loc_name, length))
	sps = sorted(list(seq.keys()))
	for sp in sps:
		tmpseq = ''.join(seq[sp])
		tmpseq = re.sub('N', '-', tmpseq)
		o.write('%s %s\n' % (groups[sp], tmpseq))
o.write(';End;\n')

#haps = {}
#for group in groups:
#	haps[group] = []
#	for ind in groups[group]:
#		haps[group].append('%s_1' % ind)
#		haps[group].append('%s_2' % ind)
#
#inds = ['%s:%s' % (group, ','.join(haps[group])) for group in haps]
#mapping = '<' + '; '.join(inds) + '>'

o.write('BEGIN PHYLONET;\n')
o.write('MCMC_SEQ -dir %s -cl 30000000 -bl 3000000 -sf 5000;\n' % (outdir))
o.write('END;\n')
o.close()
