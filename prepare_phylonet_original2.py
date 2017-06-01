import re
import glob
import os
import random
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", required=True, help="Fle to run.")
args = parser.parse_args()

num_loci = 200
minpic = 0.005
minlen = 150
miss = 0.9
numruns = 20

def get_seq(locus):
	f = open(locus, 'r')
	s = {}
	length = None
	seqid = ''
	for l in f:
		if re.search('>', l):
			seqid = re.search('>(\S+)', l).group(1)
			s[seqid] = ''
		else:
			s[seqid] += l.rstrip().upper()
	f.close()

	for ind, seq in s.items():
		length = len(seq)

	return s, length					


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


def get_pics(seq, seqlen):
	n_pics = 0
	for i in range(0, seqlen):
		bases = [seq[id][i] for id in seq]
		n_pics += pic(bases)
	picval= n_pics / float(seqlen)

	return picval


clusters = {    'coggeri': ['Lampropholis_amicula', 'Lampropholis_coggeri_C', 
				'Lampropholis_coggeri_N', 'Lampropholis_coggeri_S', 
				'Lampropholis_coggeri_EU'],
                'robertsi': ['Lampropholis_robertsi_BFHR', 'Lampropholis_robertsi_BK',
	 		     'Lampropholis_robertsi_CU', 'Lampropholis_robertsi_TU', 
			     'Lampropholis_amicula'],
                'carlia': ['Carlia_rhomboidalis_N', 'Carlia_rhomboidalis_S', 
			   'Carlia_rubrigularis_N', 'Carlia_rubrigularis_S', 
                           'Carlia_storri', 'Carlia_wandalthini']
        }

d = pd.read_csv(args.file)
groups = {}
for ix, row in d.iterrows():
	groups[row['sample']] = row['lineage']

for cluster in ['robertsi', 'carlia', 'coggeri']:
	locidir = os.path.join('/scratch/drabosky_flux/sosi/AWT_delimit/original_approach/', cluster, 'loci')
	loci = glob.glob(locidir + '/*aln')
	allseq = {}
	for ix, locus in enumerate(loci):
		if ix % 100 == 0:
			print(cluster, ix)
		locname = re.search('(ENS.*).fasta', locus).group(1)
		seq, loclen = get_seq(locus)
		# print(loclen)
		# check length
		if loclen >= minlen:
			# check missingness
			sps = []
			for ind in seq:
				newind = re.sub('_\d$', '', ind)
				newind = re.sub('_R_', '', newind)
				sps.append(groups[newind])
			sps = list(set(sps))
			# print((len(sps) / float(len(clusters[cluster]))))
			if (len(sps) / float(len(clusters[cluster]))) >= miss:
				picval = get_pics(seq, loclen)
				# print(picval)
				# check pics
				if picval >= minpic:
					# winner winner, chicken dinner!
					# initialize hash
					allseq[locname] = {}

					tmp = {}
					for ind, s in seq.items():
						newind = re.sub('_\d$', '', ind)
						newind = re.sub('_R_', '', newind)
						sp = groups[newind]
						if sp not in tmp and sp in clusters[cluster]:
							tmp[sp] = []
						tmp[sp].append(s)
					# randomly pick one haplo per locus
					for sp in tmp:
						allseq[locname][sp] = random.choice(tmp[sp])

	for ix in range(0, numruns):
		print(cluster, ix)
		sel = random.sample(list(allseq.keys()), num_loci)

		tot_len = sum([len(list(allseq[loc].values())[0]) for loc in sel])
		num_inds = len(clusters[cluster])

		outdir = os.path.join('/scratch/drabosky_flux/sosi/AWT_delimit/original_approach/', cluster, 'run%s' % ix)
		if not os.path.isdir(outdir):
			os.mkdir(outdir)

		out = os.path.join(outdir, 'phylonet_miss%s_pic%s.nex' % (miss, minpic))
		o = open(out, 'w')
		o.write('#NEXUS\n')
		o.write('Begin data;\n')
		o.write('\tDimensions ntax=%s nchar=%s;\n' % (num_inds, tot_len))
		o.write('\tFormat datatype=dna symbols="ACTG" missing=? gap=-;\n')
		o.write('\tMatrix\n')

		for loc_name in sel:
			seq = allseq[loc_name]
			length = len(list(seq.values())[0])
			o.write('[%s, %s]\n' % (loc_name, length))
			# print('[%s, %s]' % (loc_name, length))
			sps = sorted(list(seq.keys()))
			for sp in sps:
				# print(len(allseq[loc_name][sp]))
				# tmpseq = re.sub('N', '-', allseq[loc_name][sp])
				o.write('%s %s\n' % (sp, allseq[loc_name][sp]))
		o.write(';End;\n')
		o.write('BEGIN PHYLONET;\n')
		o.write('MCMC_SEQ -dir %s -cl 30000000 -bl 3000000 -sf 5000;\n' % (outdir))
		o.write('END;\n')
		o.close()
