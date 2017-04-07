import argparse
import glob
import os
import pandas as pd
import re
import subprocess
import random

"""
Sonal Singhal
created on 30 March 2017
Written assuming nothing!
"""

def get_args():
        parser = argparse.ArgumentParser(
                        description="This creates the nexus files for STARBEAST.",
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter
                        )

        # file
        parser.add_argument(
                '--file',
                type=str,
                default=None,
                help='File with information for phylogeny making.'
                )

        # output dir
        parser.add_argument(
                '--outdir',
                type=str,
                default=None,
                help='Output directory for phylogeny if not '
                     'running in context of pipeline.'
                )
	
        # run number
        parser.add_argument(
                '--run',
                type=int,
                default=None,
                help='Run number, to allow multiple tests'
                )

	return parser.parse_args()


def get_sp_loci(args):
        d = pd.read_csv(args.file)
        d = d.groupby('lineage')
        groups = {}
        for lineage, group in d:
                groups[lineage] = group['sample'].tolist()

	outdir = os.path.join(args.outdir, 'starbeast')
	if not os.path.isdir(outdir):
		os.mkdir(outdir)
	subdir = os.path.join(outdir, 'run%s' % args.run)
	if not os.path.isdir(subdir):
		os.mkdir(subdir)
	
	indir = os.path.join(args.outdir, 'partitionfinder', 'run%s' % args.run)

	return groups, subdir, indir


def get_seq(seqfile):
	
	f = open(seqfile, 'r')

	id = ''
	s = {}

	# first line of phylip file
	head = f.readline()

	for l in f:
		d = re.split('\s+', l.rstrip())
		s[d[0]] = d[1]
	f.close()

	return s


def get_partitionfinder(indir):
	infile = os.path.join(indir, 'analysis/best_scheme.txt')
	i = open(infile, 'r')

	sets = {}
	read = 0
	for l in i:
		if re.search('\#nexus', l):
			read = 1
		if read:
			if re.search('charset', l):
				name = re.search('(Subset\d+)', l.strip()).group(1)
				subsets = re.findall('(\d+\-\d+)', l.rstrip())
				sets[name] = {'coords': [], 'model': None}

				for subset in subsets:
					subset = re.findall('(\d+)', subset)
					subset = [int(x) - 1 for x in subset]
					# convert to index python
					sets[name]['coords'].append(subset)
			if re.search('charpartition', l):
				read = 0
				models = re.search('= (.*);', l.rstrip()).group(1)
				models = re.split(',\s+', models)
				for model in models:
					model = re.split(':', model)
					sets[model[1]]['model'] = model[0]

	return sets


def print_output(groups, outdir, indir, sets):
	seqfile = glob.glob(indir + '/*phy')[0]
	seq = get_seq(seqfile)

	# write out ind - species file
	outfile1 = os.path.join(outdir, 'species_inds.txt')
	o = open(outfile1, 'w')
	o.write('traits\tspecies\n')
	for group, inds in groups.items():
		for ind in inds:
			o.write('%s_1\t%s\n' % (ind, group))
			o.write('%s_2\t%s\n' % (ind, group))
	o.close()

	# write out models for ref
	outfile2 = os.path.join(outdir, 'models.txt')
	o = open(outfile2, 'w')
	models = sorted(list(sets.keys()))
	for model in models:
		o.write('%s\t%s\n' % (model, sets[model]['model']))
	o.close()

	ix = 0
	# make dir for nexus
	subdir = os.path.join(outdir, 'nexus')
	if not os.path.isdir(subdir):
		os.mkdir(subdir)
	for loci in sets:
		# creates the new concatenated alignment
		newseq = {}
		for s in seq:
			newseq[s] = ''
		for subset in sets[loci]['coords']:
			start = subset[0]
			end = subset[1] + 1
			for s in seq:
				newseq[s] += seq[s][start:end]
		
		numseq = len(newseq)
		loci_len = len(list(newseq.values())[0])
		out = os.path.join(subdir, '%s.nex' % loci)
		o = open(out, 'w')
		o.write('#NEXUS\n')
		o.write('Begin data;\n')
		o.write('Dimensions ntax=%s nchar=%s\n;' % (numseq, loci_len))
		o.write('Format datatype=dna missing=N gap=-;\n')
		o.write('Matrix\n')
		for ind, s in newseq.items():
			o.write('%s\t%s\n' % (ind, s))
		o.write(';\n')
		o.write('End;\n')
		o.write('BEGIN SETS;\n')
		cur_start = 1
		for subset in sets[loci]['coords']:
			loci_len = subset[1] - subset[0] + 1
			end = cur_start + loci_len - 1
			o.write('charset set%s = %s-%s;\n' % (ix, cur_start, end))
			cur_start = end + 1
			ix += 1
		o.write('END;\n')
		o.close()
		

def main():
	args = get_args()
	groups, outdir, indir = get_sp_loci(args)
	sets = get_partitionfinder(indir)
	print_output(groups, outdir, indir, sets)

if __name__ == "__main__":
	main()
