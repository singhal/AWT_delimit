import argparse
import glob
import os
import pandas as pd
import re
import subprocess
import random

"""
Sonal Singhal
created on 23 June 2016
Written assuming nothing!
"""

def get_args():
        parser = argparse.ArgumentParser(
                        description="This creates the files that then get " 
                                    "aligned in the next script.",
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter
                        )

        # file
        parser.add_argument(
                '--file',
                type=str,
                default=None,
                help='File with information for phylogeny making.'
                )

        # miss
        parser.add_argument(
                '--miss',
                type=float,
                default=None,
                help='How much missing data will you tolerate?'
                )

        # output dir
        parser.add_argument(
                '--outdir',
                type=str,
                default=None,
                help='Output directory for phylogeny if not '
                     'running in context of pipeline.'
                )

	# number of loci
        parser.add_argument(
                '--nloci',
                type=int,
                default=None,
                help='Number of loci to use'
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

	outdir = args.outdir

	loc_file = os.path.join(outdir, 'locus_data.csv')
	d = pd.read_csv(loc_file)

	loci = d.ix[d.missingness >= args.miss, 'locus'].tolist()
	loci = random.sample(loci, args.nloci)

	return outdir, groups, loci


def make_concatenated(args, outdir, sps, loci):
	subdir = os.path.join(outdir, 'partitionfinder')
	if not os.path.isdir(subdir):
		os.mkdir(subdir)

	subdir = os.path.join(subdir, 'run%s' % args.run)
        if not os.path.isdir(subdir):
                os.mkdir(subdir)

	# where the alignments are
	seqdir = os.path.join(outdir, 'alignments')

	file1 = os.path.join(subdir, 'concatenated%s.phy' % 
			    (args.miss))
	file2 = os.path.join(subdir, 'partition_finder.cfg')

	seq = {}
	for group, inds in sps.items():
		for ind in inds:
			seq[ind + '_1'] = ''
			seq[ind + '_2'] = ''
	partitions = {}
	cur_pos = 1
	
	for locus in loci:
		f1 = os.path.join(seqdir, '%s.fasta.aln' % locus)

		f = open(f1, 'r')
		id = ''
		s = {}
		for l in f:
			if re.search('>', l):
				id = re.search('>(\S+)', l).group(1)
				# get rid of reverse
				id = re.sub('^_R_', '', id)

				s[id] = ''
			else:
				s[id] += l.rstrip()
		f.close()

		for sp, tmpseq in s.items():
			s[sp] = re.sub('\s+', '', tmpseq)

		loc_length = len(s[s.keys()[0]])
		partitions[locus] = [cur_pos, (cur_pos + loc_length - 1)]
		cur_pos = cur_pos + loc_length

		null = '-' * loc_length
		for sp in seq:
			if sp not in s:
				seq[sp] += null
			else:
				seq[sp] += s[sp]

	f = open(file1, 'w')
	f.write(' %s %s\n' % (len(seq), len(seq[seq.keys()[0]])))
	for sp, s in seq.items():
		f.write('%s   %s\n' % (sp, s))
	f.close()

	f = open(file2, 'w')
	f.write("alignment = concatenated%s.phy;\n" % args.miss)
	f.write("branchlengths = linked;\n")
	f.write("models = all;\n")
	f.write("model_selection = aicc;\n\n")

	f.write("[data_blocks]\n")
        for locus in loci:
                f.write('%s = %s-%s;\n' % (locus, partitions[locus][0], 
				           partitions[locus][1]))
	f.write("\n")
	f.write("[schemes]\n")
	f.write("search=greedy;\n")

        f.close()


def main():
	args = get_args()
	outdir, sp, loci = get_sp_loci(args)
	make_concatenated(args, outdir, sp, loci)

if __name__ == "__main__":
	main()
