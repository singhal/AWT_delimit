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

        # coords file
        parser.add_argument(
                '--cfile',
                type=str,
                default=None,
                help='File with coordinates.'
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

        # subsample
        parser.add_argument(
                '--subsample',
                default=False,
                action='store_true',
                help='Include flag if you only want to include '
                     'one haplotype per lineage; useful for '
                     'downstream concatenated analyses'
                )


	return parser.parse_args()


def get_sp_loci(args):
        d = pd.read_csv(args.file)
        d = d.groupby('lineage')
        groups = {}
        for lineage, group in d:
                groups[lineage] = group['sample'].tolist()

	if args.subsample:
		for lineage, group in groups.items():
			groups[lineage] = [random.choice(groups[lineage])]
	outdir = args.outdir

	loc_file = os.path.join(outdir, 'locus_data.csv')

	d1 = pd.read_csv(loc_file)
	d2 = pd.read_csv(args.cfile)

	d = d1.merge(d2, left_on="locus", right_on="exon")
	d = d[d.coords.notnull()]
	d = d[d.missingness >= args.miss]	
	loci = d.locus.tolist()	
	loci = random.sample(loci, args.nloci)

	return outdir, groups, loci, d


def make_concatenated(args, outdir, sps, loci, d):
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
			if not args.subsample:
				seq[ind + '_2'] = ''
	partitions = {}
	cur_pos = 1
	
	for locus in loci:
		f1 = os.path.join(seqdir, '%s.fasta.aln' % locus)
		coords = d.ix[d.locus == locus, 'coords'].tolist()[0]
		coords = [int(x) for x in re.split('_', coords)]
		# just keep the first annotated exon
		# if there are more than one
		# minus one because of indexing for both
		# add one to end because of python handling of lists
		c = [coords[0] - 1, coords[1]] 

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
			tmpseq = list(re.sub('\s+', '', tmpseq))
			tmpseq = tmpseq[c[0]:c[1]]
			s[sp] = ''.join(tmpseq)

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
                f.write('%s_1 = %s-%s\\3;\n' % (locus, partitions[locus][0], 
				           partitions[locus][1]))
		f.write('%s_2 = %s-%s\\3;\n' % (locus, partitions[locus][0] + 1,
                                           partitions[locus][1]))
		f.write('%s_3 = %s-%s\\3;\n' % (locus, partitions[locus][0] + 2,
                                           partitions[locus][1]))
	f.write("\n")
	f.write("[schemes]\n")
	if args.nloci > 200:
		f.write("search=rcluster;\n")
	else:
		f.write("search=greedy;\n")

        f.close()


def main():
	args = get_args()
	outdir, sp, loci, d = get_sp_loci(args)
	make_concatenated(args, outdir, sp, loci, d)

if __name__ == "__main__":
	main()
