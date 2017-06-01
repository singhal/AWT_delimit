import re
import pandas as pd
import glob
import os
import multiprocessing as mp
import subprocess

clusters = {	'coggeri': ['Lampropholis_amicula', 'Lampropholis_coggeri_C', 'Lampropholis_coggeri_N', 'Lampropholis_coggeri_S', 'Lampropholis_coggeri_EU'],
		'robertsi': ['Lampropholis_robertsi_BFHR', 'Lampropholis_robertsi_BK', 'Lampropholis_robertsi_CU', 'Lampropholis_robertsi_TU', 'Lampropholis_amicula'],
		'carlia': ['Carlia_rhomboidalis_N', 'Carlia_rhomboidalis_S', 'Carlia_rubrigularis_N', 'Carlia_rubrigularis_S', 'Carlia_storri', 'Carlia_wandalthini']
	}

'''
exons = {}
f = open('/scratch/drabosky_flux/sosi/AWT_delimit/targetExons.fa', 'r')
for l in f:
	if re.search('>', l):
		exon = re.search('>(\S+)', l).group(1)
		exon = re.sub('_[A-Z|a-z]+$', '', exon)
		exons[exon] = 1
f.close()

def get_seq(seqfile, seq):
        o = open(seqfile, 'r')
        seqid = ''
        for l in o:
                if re.search('>', l):
                        seqid = re.search('>(\S+)', l).group(1)
                        seq[seqid] = ''
                else:
                        seq[seqid] += l.rstrip()

        o.close()

        return seq

for cluster in clusters:
	outdir = os.path.join('/scratch/drabosky_flux/sosi/AWT_delimit/original_approach', cluster, 'loci')
	if not os.path.isdir(outdir):
		os.mkdir(outdir)

	for exon in exons:
		seq = {}
		for sp in clusters[cluster]:
			seqfile = os.path.join('/scratch/drabosky_flux/sosi/AWT_delimit/original_approach', 'variants', sp, '%s.aln.fasta' % exon)
			if os.path.isfile(seqfile):
				seq = get_seq(seqfile, seq)

		keep = {}
		for ind, s in seq.items():
			s = re.sub('n', '', s)
			s = re.sub('N', '', s)
			if len(s) >= 50:
				keep[ind] = s
		
		if len(keep) > 3:
			out = os.path.join('/scratch/drabosky_flux/sosi/AWT_delimit/original_approach', cluster, 'loci', '%s.fasta' % exon)
			o = open(out, 'w')
			for ind, s in keep.items():
				o.write('>%s\n%s\n' % (ind, s))
			o.close()
'''

def align(file):
      
        aln_out = file.replace('.fasta', '.fasta.aln')
        proc = subprocess.call("mafft --maxiterate 1000 --globalpair "
                               "--adjustdirection --quiet %s > %s" %
                               (file, aln_out), shell=True)

        # os.remove(file)
        return aln_out


mafft = 'mafft'
for cluster in clusters:

	cdir = os.path.join('/scratch/drabosky_flux/sosi/AWT_delimit/original_approach', cluster, 'loci')
	files = glob.glob(cdir + '/*fasta')

	pool = mp.Pool(16)
	alns = pool.map(align, files)
