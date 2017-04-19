import re
import glob
import pandas as pd
import numpy as np
import subprocess
import os

# nucleotides to proteins
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

# target file
ann = '/scratch/drabosky_flux/sosi/AWT_delimit/targetExons.fa'

def translate(seq, gencode):
	seq = [gencode[''.join(seq[i:i+3])] 
			if ''.join(seq[i:i+3]) in gencode else 'X' 
			for i in range(0, int(len(seq)/3)*3, 3)]
	seq = ''.join(seq)
	
	return(seq)


def get_seq(seqfile):
	seq = {}
	id = ''
	
	f = open(seqfile, 'r')
	for l in f:
		l = l.rstrip()
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			seq[id] = ''
		else:
			seq[id] += l.upper()
	f.close() 

	return seq

def get_frame(ann, gencode):
	drop = []
	allframes = {}

	for id, s in ann.items():
		s = list(s)
		frames = {}
		for frame in [0, 1, 2]:
			cds_len = len(s) - 1
			num_aa = int(cds_len / 3.0)
	
			cds = s[frame:(1 + frame + 3 * num_aa)]
			aa = translate(cds, gencode)
			if aa.count('*') == 0:
				frames[frame] = 1

		id1 = re.sub('_[A-Z|a-z]+$', '', id)

		if len(frames) != 1:
			if id1 not in drop:
				drop.append(id1)
		else:
			allframes[id1] = [id, frames.keys()[0]]

	return allframes, drop

target = get_seq(ann)
allframes, drop = get_frame(target, gencode)

align = '/scratch/drabosky_flux/sosi/AWT_delimit/phylogeny/coordinates.csv'
a = open(align, 'w')
a.write('exon,coords\n')

seqfiles = glob.glob('/scratch/drabosky_flux/sosi/AWT_delimit/phylogeny/alignments/*aln')
for seqfile in seqfiles:
	ann = get_seq(seqfile)
	locname = re.sub('.*alignments/', '', seqfile)
	locname = re.sub('.fasta.aln', '', locname)

	if locname not in drop:
		aaseq = target[allframes[locname][0]]
		aaseq = list(aaseq)
		aa = translate(aaseq[allframes[locname][1]:], gencode)
		
		p_out = '%s_prot.fasta' % (locname)
		o = open(p_out, 'w')
		o.write('>prot\n%s\n' % aa)
		o.close()
		
		d_out = '%s_dna.fasta' % (locname)
		o = open(d_out, 'w')
		seq = re.sub('-', 'N', ann.values()[0])
		o.write('>dna\n%s\n' % seq)
		o.close()
		
		call = subprocess.Popen("exonerate -m protein2genome -q %s -t %s --showtargetgff TRUE --showalignment NO --showvulgar 0 -n 1 | grep 'cds'" % (p_out, d_out), shell=True, stdout=subprocess.PIPE)
		lines = [line.rstrip() for line in call.stdout]
		coords = []
		for line in lines:
			d = re.split('\s+', line)
			coords.append(int(d[3]))
			coords.append(int(d[4]))
		coords = list(set(coords))
		coords = [str(x) for x in sorted(coords)]
		os.remove(d_out)
		os.remove(p_out)
		if len(coords) > 1:
			a.write('%s,%s\n' % (locname, '_'.join(coords)))
		else:
			a.write('%s,NA\n' % (locname))
	else:
		a.write('%s,NA\n' % (locname))
a.close()



