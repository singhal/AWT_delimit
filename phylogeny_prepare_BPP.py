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

        # species tree
        parser.add_argument(
                '--sp',
                type=str,
                default=None,
                help='Species tree'
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

def define_ctl(subdir, file1, file2, sps, sptree, nloci, theta1, theta2, tau1, tau2):
	outfile1 = os.path.join(subdir, 'out_theta%s_%s_tau%s_%s.txt' % 
			(theta1, theta2, tau1, tau2))
	outfile2 = os.path.join(subdir, 'out_theta%s_%s_tau%s_%s_mcmc.txt' % 
                        (theta1, theta2, tau1, tau2))
	outfile = os.path.join(subdir, 'theta%s_%s_tau%s_%s.ctl.txt' %
                        (theta1, theta2, tau1, tau2))

	random_num = random.randint(-100, 100)
	n_sp = len(sps)
	sp_names = ' '.join(sorted(sps.keys()))
	n_hap = [2 * len(sps[sp]) for sp in sps]
	n_hap = ' '.join([str(x) for x in n_hap])

	info = [
  	"seed = %s" % random_num, 
  	"seqfile = %s" % file1,
  	"Imapfile = %s" % file2,
  	"outfile = %s" % outfile1,
  	"mcmcfile = %s" % outfile2, 
  	"speciesdelimitation = 1 0 1",
  	"speciestree = 0",
  	"species&tree = %s %s" % (n_sp, sp_names),
  	                "   %s" % (n_hap),
  	                "   %s" % sptree,
  	"usedata = 1",
  	"nloci = %s" % nloci,
  	"cleandata = 0",
  	"thetaprior = %s %s" % (theta1, theta2),
  	"tauprior = %s %s" % (tau1, tau2),
  	"locusrate = 1 2.0",
  	"finetune = 1: 0.5 0.0015 0.0006 0.0004 0.06 0.2 1.0",
  	"print = 1 0 0 0",
  	"burnin = 10000",
  	"sampfreq = 5",
  	"nsample = 100000"
  	]

	return info, outfile

def make_concatenated(args, outdir, sps, loci):
	subdir = os.path.join(outdir, 'BPP')
	if not os.path.isdir(subdir):
		os.mkdir(subdir)

	subdir = os.path.join(subdir, 'run%s' % args.run)
        if not os.path.isdir(subdir):
                os.mkdir(subdir)

	# where the alignments are
	seqdir = os.path.join(outdir, 'alignments')

	file1 = os.path.join(subdir, 'loci_nloci%s_miss%s.phy' % 
			    (args.nloci, args.miss))
	file2 = os.path.join(subdir, 'species_imap.txt')
	f = open(file2, 'w')
	groups = {}
	for group, inds in sps.items():
		for ind in inds:
			groups[ind] = group
			f.write('%s_1 %s\n'  % (ind, group))
			f.write('%s_2 %s\n'  % (ind, group))
	f.close()
	
	o = open(file1, 'w')
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

		numseq = len(s)
		loc_len = len(list(s.values())[0])

		o.write('     %s     %s\n' % (numseq, loc_len)) 
		for id, seq in s.items():
			o.write('%s^%s  %s\n' % (id, id, seq))
		o.write('\n\n')


	theta = [(2, 2000), (1, 10), (1, 10)]
	tau = [(2, 2000), (1, 10), (2, 2000)]
	for th, ta in zip(theta, tau):
		info, out = define_ctl(subdir, file1, file2, sps,
					args.sp, args.nloci, 
					th[0], th[1], ta[0], ta[1])
		o = open(out, 'w')
		for i in info:
			o.write(i + '\n')
		o.close()


def main():
	args = get_args()
	outdir, sp, loci = get_sp_loci(args)
	make_concatenated(args, outdir, sp, loci)

if __name__ == "__main__":
	main()
