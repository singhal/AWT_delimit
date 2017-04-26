import re
import argparse
import pandas as pd
import os
import numpy as np
import gzip
import pickle

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", required=True, help='File with the sample data.')
parser.add_argument("-n", "--name", required=True, help='Output name.')
args = parser.parse_args()

indir = '/Users/Sonal/Dropbox/Carlia_Lampropholis_Species_Delimitation/data/'
outfile = os.path.join(indir, '%s_divergence.csv' % args.name)
out = open(outfile, 'w')
out.write('sp1,sp2,metric,denom,value\n')

def site_count(gencode):

	sites = {}
	for codon in gencode:
		codon_sites = []
		for ix, base in enumerate(codon):
			tmp = list(codon)
			aa = []
			for mut in ['A', 'T', 'C', 'G']:
				tmp[ix] = mut
				aa.append(gencode[''.join(tmp)])
			num = len(set(aa))
			codon_sites.append( ( 4 - num) / 3.0)
		sites[codon] = codon_sites

	return sites

def fst_reich(counts, sample_sizes):
	counts1 = np.array(counts)[0]
	counts2 = np.array(counts)[1]

	sample_sizes1 = np.array(sample_sizes).astype('float')[0]
	sample_sizes2 = np.array(sample_sizes).astype('float')[1]

	h1 = counts1 * (sample_sizes1 - counts1) / (sample_sizes1 * (sample_sizes1 - 1))
	h2 = counts2 * (sample_sizes2 - counts2) / (sample_sizes2 * (sample_sizes2 - 1))
	
	N = []
	D = []

	for _a1, _a2, _n1, _n2, _h1, _h2 in zip(counts1, counts2, sample_sizes1, sample_sizes2, h1, h2):
		n = ((_a1 / _n1) - (_a2 / _n2)) ** 2 - (_h1 / _n1) - (_h2 / _n2)
		N.append(n)
		d = n + _h1 + _h2
		D.append(d)

	F = np.sum(N) / np.sum(D)

	return F

def get_groups(sfile):
	d = pd.read_csv(sfile)

	sps = {}
	inds = {}

	for ix, row in d.iterrows():
		ind = row['sample']
		sp = row['lineage']

		if sp not in sps:
			sps[sp] = []

		sps[sp].append(ind)
		inds[ind] = sp

	for sp, ind in sps.items():
		sps[sp] = sorted(ind)

	return sps, inds

def get_variants(indir, sps):
	var = {}

	for sp in sps:
		vcf = os.path.join(indir, 'variants', 
			'%s.qual_filtered.cov_filtered.vcf.gz' % sp)
		f = gzip.open(vcf, 'r')

		for l in f:
			l = l.decode('utf-8')
			if re.search('CHROM', l):
				ids = re.split('\t', l.rstrip())[9:]
			elif not re.search('#', l) and not re.search('INDEL', l):
				d = re.split('\t', l.rstrip())

				c = d[0]
				pos = int(d[1])
				alleles = [d[3]] + re.split(',', d[4])

				if c not in var:
					var[c] = {}
				if pos not in var[c]:
					var[c][pos] = {}
				
				snpcode = {}
				for ix, a in enumerate(alleles):
					snpcode[str(ix)] = a
				snpcode['.'] = 'N'

				genos = d[9:]
				genos = [re.search('^(\S/\S)', x).group(1) \
							for x in genos]
				for id, geno in zip(ids, genos):
					geno = re.split('/', geno)
					snps = [snpcode[g] for g in geno]
					var[c][pos][id] = snps

	return var

def get_coords(indir, sp):
	cfile = os.path.join(indir, 'PRG', '%s_coordinates.csv' % sp)
	c = {}

	f = open(cfile, 'r')
	head = f.readline()

	# only do the first coding sequence 
	# to make life easier
	for l in f:
		d = re.split(',', l)
		if d[1] != 'NA':
			pos = re.findall('(\d+)', d[1])
			
			c[d[0]] = [int(pos[0]), int(pos[1])]

	return c

def calc_pi_prop(alleles):
	alleles = dict([(x, alleles.count(x)) for x in set(alleles)])
	if len(alleles) > 1:
		# tot values
		n = float(np.sum(list(alleles.values())))
		# minor count
		j = float(np.min(list(alleles.values())))
		pi_prop = (2 * j * (n - j)) / (n * (n - 1))
	else:
		pi_prop = 0

	return pi_prop

def calc_pi(sp, sps, var, c):
	ids = sps[sp]
	denom = 0
	diff = 0

	for exon in var:
		for pos in var[exon]:
			# within coding sequence
			if pos >= c[exon][0] and pos <= c[exon][1]:
				alleles = []
				for id in ids:
					if id in var[exon][pos]:
						alleles += var[exon][pos][id]
				alleles = [a for a in alleles if a != 'N']

				# need at least two chromosomes
				if len(alleles) > 1:
					
					alleles = dict([(x, alleles.count(x)) for x in set(alleles)])
					# don't mess with multiallelics
					if len(alleles) < 3:
						denom += 1
						diff += calc_pi_prop(alleles)

	pi = np.nan
	if denom > 0:
		pi = diff / float(denom)

	# two times length of inds because diploids!
	return pi, denom, 2 * len(ids)

def get_alleles(ids, var2):
	alleles = []
	for id in ids:
		if id in var2:
			alleles += var2[id]
	alleles = [a for a in alleles if a != 'N']
	return(alleles)

def calc_div_prop(a1, a2):
	n_c = 0
	n_diff = 0
	for x in a1:
		for y in a2:
			n_c += 1
			if x != y:
				n_diff += 1

	return n_diff / float(n_c)

def calc_div_sub(sp1, sp2, sps, var, c):
	ids1 = sps[sp1]
	ids2 = sps[sp2]

	denom = 0
	diff = 0
	raw_pi1 = 0
	raw_pi2 = 0

	for exon in var:
		for pos in var[exon]:
			# within coding sequence
			if pos >= c[exon][0] and pos <= c[exon][1]:
				a1 = get_alleles(ids1, var[exon][pos])
				a2 = get_alleles(ids2, var[exon][pos])

				# need at least two chromosomes
				if len(a1) > 1 and len(a2) > 1:
					# don't mess with multiallelics
					if len(set(a1)) < 3 and len(set(a2)) < 3:
						denom += 1

						diff += calc_div_prop(a1, a2)
						raw_pi1 += calc_pi_prop(a1)
						raw_pi2 += calc_pi_prop(a2)

	pi_btn = pi1 = pi2 = np.nan
	if denom > 0:
		pi_btn = diff / float(denom)
		pi1 = raw_pi1 / float(denom)
		pi2 = raw_pi2 / float(denom)

	return pi_btn, pi1, pi2, denom

def calc_div(sps, var, c):
	species = sorted(list(sps.keys()))

	for ix, sp1 in enumerate(species):
		for sp2 in species[(ix+1):]:

			pi_btn, pi1, pi2, denom = calc_div_sub(sp1, sp2, sps, var, c)
			pi_net = pi_btn - (pi1 + pi2) * 0.5

			out.write('%s,NA,pi,%s,%.5f\n' % (sp1, denom, pi1))
			out.write('%s,NA,pi,%s,%.5f\n' % (sp2, denom, pi2))
			out.write('%s,%s,d_xy,%s,%.5f\n' % (sp1, sp2, denom, pi_btn))
			out.write('%s,%s,d_a,%s,%.5f\n' % (sp1, sp2, denom, pi_net))
			
def calc_fst_sub(sp1, sp2, sps, var, c):
	ids1 = sps[sp1]
	ids2 = sps[sp2]

	counts = {sp1: [], sp2: []}
	sizes = {sp1: [], sp2: []}

	for exon in var:
		for pos in var[exon]:
			# within coding sequence
			if pos >= c[exon][0] and pos <= c[exon][1]:
				a1 = get_alleles(ids1, var[exon][pos])
				a2 = get_alleles(ids2, var[exon][pos])

				# only want to work with positions where 
				# there are two snps and two snps only
				if len(set(a1 + a2)) == 2:
					# this is the variable site
					# determined arbitarly
					# fst doesn't have minor / major ?? (i think?)
					allele = list(set(a1 + a2))[0] 

					if len(a1) == 0:
						count1 = np.nan
					else:
						count1 = a1.count(allele)
					counts[sp1].append(count1)
					sizes[sp1].append(len(a1))

					if len(a2) == 0:
						count2 = np.nan
					else:
						count2 = a2.count(allele)
					counts[sp2].append(count2)
					sizes[sp2].append(len(a2))

	# confirm sample sizes
	alleles = np.array([counts[sp1], counts[sp2]])
	sizes = np.array([sizes[sp1], sizes[sp2]])
	to_mask = np.any(np.isnan(alleles), axis=0)
	alleles = alleles[:, -to_mask]
	sizes = sizes[:, -to_mask]

	if len(alleles[0]) > 0:
		fst = fst_reich(alleles, sizes)
	else:
		fst = np.nan

	return fst, len(alleles[0])

def calc_fst(sps, var, c):
	species = sorted(list(sps.keys()))

	for ix, sp1 in enumerate(species):
		for sp2 in species[(ix+1):]:

			fst, denom = calc_fst_sub(sp1, sp2, sps, var, c)
			out.write('%s,%s,F_ST,%s,%.5f\n' % (sp1, sp2, denom, fst))

def get_seq(indir, sp):
	seqfile = os.path.join(indir, 'PRG', '%s.fasta' % sp)
	f = open(seqfile, 'r')

	seq = {}
	id = ''
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			seq[id] = ''
		else:
			seq[id] += l.rstrip()

	return seq

def calc_div_silent(sps, var, coords, gencode, sites):
	
	species = sorted(list(sps.keys()))

	for ix, sp1 in enumerate(species):
		ids1 = sps[sp1]
		seq = get_seq(indir, sp1)
		for sp2 in species[(ix+1):]:
			res = {'denom': 0, 'pi1': 0, 'pi2': 0, 'btn': 0}
			ids2 = sps[sp2]
			for c in coords:
				start = coords[c][0]
				end = coords[c][1]
				cds = seq[c][(start - 1):end]

				for i in range(start, end, 3):
					rel_pos = i - start
					codon = cds[rel_pos: rel_pos + 3]
					pos = [i, i+1, i+2]

					if codon in gencode:
						vals = sites[codon]
						for ix, p in enumerate(pos):
							if p in var[c]:

								alleles1 = get_alleles(ids1, var[c][p])
								alleles2 = get_alleles(ids2, var[c][p])

								# make sure at least two chrs in sp1
								# make sure at least two chrs in sp2
								# make sure just one or two alleles 
								if len(alleles1) > 1 and len(alleles2) > 1:
									alleles = list(set(alleles1 + alleles2))
									if len(alleles) < 3:
										res['denom'] += vals[ix]

										if len(alleles) == 2:
											syn = True
											mut1 = list(codon)
											mut2 = list(codon)

											mut1[ix] = alleles[0]
											mut2[ix] = alleles[1]

											aa1 = gencode[''.join(mut1)]
											aa2 = gencode[''.join(mut2)]

											if aa1 != aa2:
												syn = False
								
											if syn:
												res['pi1'] += calc_pi_prop(alleles1)
												res['pi2'] += calc_pi_prop(alleles2)
												res['btn'] += calc_div_prop(alleles1, alleles2)
			pi1 = res['pi1'] / res['denom']
			pi2 = res['pi2'] / res['denom']
			pi_btn = res['btn'] / res['denom']
			pi_net = pi_btn - 0.5 * (pi1 + pi2)

			out.write('%s,NA,pi_silent,%s,%.5f\n' % (sp1, res['denom'], pi1))
			out.write('%s,NA,pi_silent,%s,%.5f\n' % (sp2, res['denom'], pi2))
			out.write('%s,%s,d_xy_silent,%s,%.5f\n' % (sp1, sp2, res['denom'], pi_btn))
			out.write('%s,%s,d_a_silent,%s,%.5f\n' % (sp1, sp2, res['denom'], pi_net))

gencode = {
	'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
	'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
	'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
	'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
	'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
	'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
	'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
	'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
	'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
	'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
	'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
	'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
	'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
	'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
	'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
	'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
	}
sites =  site_count(gencode)

# 0.005 to 0.02 at silent 
# get groups
sps, inds = get_groups(args.file)
# get variant data
# what should my data structure be?
# exon, pos, inds, variants?
# var = get_variants(indir, sps)
# pickle.dump(var, open( "var.pickle", "wb" ))
var = pickle.load(open( "var.pickle", "rb" ))
# get cds coordinates
# because of logic of ref approach
# coords will be same across all lineages
c = get_coords(indir, list(sps.keys())[0])
# calculate divergence btn & within species
# only within coding sequence though
calc_div(sps, var, c)
# calculate Fst
calc_fst(sps, var, c)
# calculate silent
calc_div_silent(sps, var, c, gencode, sites)

out.close()
