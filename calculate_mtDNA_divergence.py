import re
import argparse
import pandas as pd
import glob
import os
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--contact", required=True, help='Contact.')
args = parser.parse_args()

# samples
d = pd.read_csv("/scratch/drabosky_flux/sosi/AWT_delimit/mtDNA/mtDNA_samples.csv")

# dir
indir = os.path.join('/scratch/drabosky_flux/sosi/AWT_delimit/mtDNA', args.contact)

# alignment
aln = glob.glob(indir + '/*phy')[0]

# partitions
part = glob.glob(indir + '/*partitions')[0]

# outfile
outfile = os.path.join('/scratch/drabosky_flux/sosi/AWT_delimit/mtDNA', '%s_mtDNA_divergence.csv' % args.contact)

def get_seq(aln):
        f = open(aln, 'r')
        s = {}
       
        for l in f:
                if not re.search('\s+\d+', l):
                        (id, seq) = re.split('\s+', l.rstrip())
                        s[id] = seq.upper()
        f.close()

        return s


def calculate_div(seq, sp1, sp2, groups):
	inds1 = [ind for ind in groups if groups[ind] == sp1]
	inds2 = [ind for ind in groups if groups[ind] == sp2]

	div = []
	legal = ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g']

	for ind1 in inds1:
		for ind2 in inds2:
			if ind1 != ind2:
				seq1 = seq[ind1]
				seq2 = seq[ind2]

				denom = 0
				diff = 0
	
				for a, b in zip(seq1, seq2):
					if a in legal and b in legal:
						denom += 1
						if a != b:
							diff += 1
				tmpdiv = np.nan
				if denom > 0:
					tmpdiv = diff / float(denom)
				div.append(tmpdiv)

	ncomp = len(div)
	div = np.mean(div)

	return div, ncomp


def calculate_div_silent(sites, seq, sp1, sp2, groups):
        inds1 = [ind for ind in groups if groups[ind] == sp1]
        inds2 = [ind for ind in groups if groups[ind] == sp2]

	div = []

        for ind1 in inds1:
                for ind2 in inds2:
                        if ind1 != ind2:
                                seq1 = seq[ind1]
                                seq2 = seq[ind2]

				denom = 0
				diff = 0

				for i in range(0, len(seq1), 3):
					codon1 = seq1[i:(i+3)]
					codon2 = seq2[i:(i+3)]

					if codon1 in gencode and codon2 in gencode:
						denom += sum(sites[codon1])
						for ix, (a, b) in enumerate(zip(codon1, codon2)):
							if a != b:
								new = list(codon1)
								new[ix] = b
					
								if gencode[codon1] == gencode[''.join(new)]:
									diff += 1
				tmpdiv = np.nan
				if denom > 0:
					tmpdiv = diff / float(denom)
				div.append(tmpdiv)
	
	n = len(div)
	div = np.mean(div)
	
	return div, n											

def get_cds(seq, part):
	p = open(part, 'r')
	start = 1e6
	end = 0

	for l in p:
		if re.search('3[,|\n]', l):
			s = re.search('(\d+)-\d+\S*3', l).group(1)
			e = re.search('\d+-(\d+)\S*3', l).group(1)
		
			if int(s) < start:
				start = int(s)
			if int(e) > end:
				end = int(e)

	cds = {}
	for id, s in seq.items():
		cds[id] = s[(start - 1):end]

	print(start,end)

	return cds


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

o = open(outfile, 'w')
o.write('lineage1,lineage2,type,n,value\n')
# get basic info
seq = get_seq(aln)
groups = {}
for id in seq:
	sp = d.ix[d['sample'] == id, 'lineage'].tolist()[0]
	groups[id] = sp
sps = list(set(groups.values()))

# calculate div
for sp1 in sps:
	for sp2 in sps:
		div, n = calculate_div(seq, sp1, sp2, groups)
		o.write('%s,%s,all,%s,%.5f\n' % (sp1, sp2, n, div))

# calculate div - silent
cds = get_cds(seq, part)
sites =  site_count(gencode)

for sp1 in sps:
	for sp2 in sps:
		div, n = calculate_div_silent(sites, cds, sp1, sp2, groups)
		o.write('%s,%s,silent,%s,%.5f\n' % (sp1, sp2, n, div))
o.close()
