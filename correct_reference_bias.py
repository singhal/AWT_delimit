import argparse
import os
import pandas as pd
import re
import subprocess

"""
Sonal Singhal
created on 21 March 2017
Written assuming:
        * samtools 1.3.1
        * GATK 3.6
        * picard 2.4.1
        * bwa 0.7.12
"""

def get_args():
	parser = argparse.ArgumentParser(
		description="Generate reference genome based on another lineage's assembly"
                            " and lineage of interest reads, written assuming "
                            " samtools 1.3.1, GATK 3.6, picard 2.4.1, and"
                            " bwa 0.7.12",
        	formatter_class=argparse.ArgumentDefaultsHelpFormatter
		)

	# sample
	parser.add_argument(
                '--lineage',
                type=str,
                default=None,
                help='Lineage for which to run script.'
                )

	# file
	parser.add_argument(
                '--file',
                type=str,
                default=None,
                help='File with sample info.'
                )

	# basedir
	parser.add_argument(
                '--dir',
                type=str,
                default=None,
                help="Full path to base dir with reads & assemblies."
                )
        
        # bwa
	parser.add_argument(
                '--bwa',
                type=str,
                default=None,
                help='bwa executable, full path.'
                )

	# samtools
        parser.add_argument(
                '--samtools',
                type=str,
                default=None,
                help='samtools executable, full path.'
                )

	# GATK
        parser.add_argument(
                '--gatk',
                type=str,
                default=None,
                help='GATK executable, full path.'
                )

	# picard
        parser.add_argument(
                '--picard',
                type=str,
                default=None,
                help='picard executable, full path.'
                )
	
	# CPUs
	parser.add_argument(
                '--CPU',
                type=int,
                default=1,
                help='# of CPUs to use in alignment.'
               )

	# memory
        parser.add_argument(
                '--mem',
                type=int,
                default=1,
                help='Memory available, as an int, in terms of Gb.'
               )
               
	# PRG
	parser.add_argument(
                '--prg',
                type=str,
                default=None,
                help="Path to PRG to use to build new ref "
                     "genome."
                )

	return parser.parse_args()


def get_info(args):
	# gets the samples associated with the lineage
	d = pd.read_csv(args.file)
	lineage = args.lineage
	samples = d.ix[d['lineage'] == lineage, 'sample'].tolist()
	samples = sorted(samples)

	# gets reads
	reads = {}
	for sample in samples:
	        read1 = os.path.join(args.dir, 'trim_reads', '%s_R1.final.fq.gz' % sample)
	        read2 = os.path.join(args.dir, 'trim_reads', '%s_R2.final.fq.gz' % sample)
	        un = os.path.join(args.dir, 'trim_reads', '%s_unpaired.final.fq.gz' % sample)
	        reads[sample] = [read1, read2, un]

	# defines the outdir
	outdir = os.path.join(args.dir, 'reference_bias')

	# makes the outdir
	if not os.path.isdir(outdir):
		os.mkdir(outdir)

	# gets the prg to mutate into the new ref
        genome = args.prg	

	return reads, lineage, genome, outdir


def prepare_seq(args, genome):
	# does all the prep necessary for the PRG
	if not os.path.isfile(genome + '.bwt'):
		subprocess.call("%s index %s" % (args.bwa, genome), shell=True)
	if not os.path.isfile(genome + '.fai'):
		subprocess.call("%s faidx %s" % (args.samtools, genome), shell=True)
	out = re.sub('.fa.*$', '.dict', genome)
	if not os.path.isfile(out):
		subprocess.call("java -jar %s CreateSequenceDictionary R=%s O=%s" % 
                                (args.picard, genome, out), shell=True)


def align_seq(args, sample, r, seq, dir, num):
	out1a = os.path.join(dir, '%s_1.sam' % sample)
	out1as = os.path.join(dir, '%s_1.bam' % sample)
	out1b = os.path.join(dir, '%s_2.sam' % sample)
	out2a = os.path.join(dir, '%s.mateFixed.bam' % sample)
	out2b = os.path.join(dir, '%s.bam' % sample)
	out3a = os.path.join(dir, '%s_1.mateFixed.sorted.bam' %  sample)
	out3b = os.path.join(dir, '%s_2.mateFixed.sorted.bam' % sample)
	out4a = os.path.join(dir, '%s_1.rg.mateFixed.sorted.bam' % sample)
        out4b = os.path.join(dir, '%s_2.rg.mateFixed.sorted.bam' % sample)
	out5a = os.path.join(dir, '%s_1.dup.rg.mateFixed.sorted.bam' % sample)
        out5b = os.path.join(dir, '%s_2.dup.rg.mateFixed.sorted.bam' % sample)
	out6 = os.path.join(dir, '%s.dup.rg.mateFixed.sorted.bam' % sample)
	intervals = os.path.join(dir, '%s.intervals' % sample)
	out7 = os.path.join(dir, '%s.realigned.dup.rg.mateFixed.sorted%s.bam' % (sample, num))
	m_1 = os.path.join(dir, '%s_1.metrics' % sample)
	m_2 = os.path.join(dir, '%s_2.metrics' % sample)

	# need a tmpdir for when sorting BAM files
	tmpdir = os.path.join(dir, sample)
	if not os.path.isdir(tmpdir):
		os.mkdir(tmpdir)

	# align
	subprocess.call("%s mem -t %s %s %s %s > %s" % (args.bwa, args.CPU, seq, r[0], r[1], out1a), shell=True)
	subprocess.call("%s mem -t %s %s %s > %s" % (args.bwa, args.CPU, seq, r[2], out1b), shell=True)
	# fixmate
	# note that had used samtools, appears to not work properly
	subprocess.call("%s view -b %s > %s" % (args.samtools, out1a, out1as), shell=True)
	subprocess.call("java -jar %s FixMateInformation I=%s O=%s" % (args.picard, out1as, out2a), shell=True)
	subprocess.call("%s view -b %s > %s" % (args.samtools, out1b, out2b), shell=True)
	# sorted
	subprocess.call("%s sort -O bam -o %s -T %s %s" % (args.samtools, out3a, tmpdir, out2a), shell=True)
	subprocess.call("%s sort -O bam -o %s -T %s %s" % (args.samtools, out3b, tmpdir, out1b), shell=True)
	# readgroup
	subprocess.call("java -jar %s AddOrReplaceReadGroups INPUT=%s OUTPUT=%s RGLB=%s RGPL=Illumina RGPU=%s RGSM=%s" % 
                          (args.picard, out3a, out4a, sample, sample, sample), shell=True)
        subprocess.call("java -jar %s AddOrReplaceReadGroups INPUT=%s OUTPUT=%s RGLB=%s RGPL=Illumina RGPU=%s RGSM=%s" % 
                          (args.picard, out3b, out4b, sample, sample, sample), shell=True)
	# mark read duplicates
	subprocess.call("java -jar %s MarkDuplicates I=%s O=%s ASO=coordinate METRICS_FILE=%s" % 
                       (args.picard, out4a, out5a, m_1), shell=True)
	subprocess.call("java -jar %s MarkDuplicates I=%s O=%s ASO=coordinate METRICS_FILE=%s" % 
                       (args.picard, out4b, out5b, m_2), shell=True)
	subprocess.call("%s merge %s %s %s" % (args.samtools, out6, out5a, out5b), shell=True)
	# indel target
	subprocess.call("%s index %s" % (args.samtools, out6), shell=True)
	subprocess.call("java -Xmx%sg -jar %s -T RealignerTargetCreator -R %s -I %s -o %s -nt %s" % 
                        (args.mem, args.gatk, seq, out6, intervals, args.CPU), shell=True)
	# indel realigner
	subprocess.call("java -Xmx%sg -jar %s -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s" % 
			(args.mem, args.gatk, seq, out6, intervals, out7), shell=True)
	

	# remove the files
	[os.remove(x) for x in [out1a, out1as, out1b, out2a, out2b, out3a, 
                                out3b, out4a, out4b, out5a, out5b, out6, 
			        out6 + '.bai', intervals]]
	
	# remove the dir
	os.rmdir(tmpdir)

	return out7


def get_seq(seqfile):
	id = ''
	seq = {}
	
	f = open(seqfile, 'r')
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			seq[id] = ''
		else:
			seq[id] += l.rstrip()

	f.close()

	for id, s in seq.items():
		seq[id] = list(s)

	return seq


def call_snps(args, files, genome, dir, num):
        raw_vcf = os.path.join(dir, '%s.raw%s.vcf' % (args.lineage, num))

	bam = '-I ' + ' -I '.join(files)

	# makes the raw VCFs, only outputting variant SNPs
	subprocess.call("java -Xmx%sg -jar %s -T UnifiedGenotyper -R %s %s -o %s "
                        "--output_mode EMIT_VARIANTS_ONLY -nt %s"
                        % (args.mem, args.gatk, genome, bam, raw_vcf, args.CPU), shell=True)

	# get in seq
	seq = get_seq(genome)

	# mutate seq
	i = open(raw_vcf, 'r')
	for l in i:
		if not re.search('#', l):
			d = re.split('\t', l.rstrip())
			# check biallelic
			if d[3] in ['A', 'T', 'C', 'G'] and d[4] in ['A','T', 'C', 'G']:
				# get af
				genos = d[9:]
				genos = [re.search('(\S/\S)', x).group(1) for x in genos]
				genos = [re.split('/', x) for x in genos]
				genos = [x for ind in genos for x in ind]
				genos = [x for x in genos if x != '.']

				if len(genos) > 0:
					af = genos.count('1') / float(len(genos))
					# mutate seq
					if af > 0.5:
						pos = int(d[1]) - 1
						seq[d[0]][pos] = d[4]
	i.close()

	out = os.path.join(dir, '%s.%s.fasta' % (args.lineage, num))
	o = open(out, 'w')
	for id, s in seq.items():
		o.write('>%s\n%s\n' % (id, ''.join(s)))
	o.close()

	return out


def main():
	# get arguments
	args = get_args()
	reads, lineage, genome, outdir = get_info(args)
	
	# round 1
	# prep sequence
	prepare_seq(args, genome)
	# do the alignments for the first time
	bamfiles1 = []
	for sample, r in reads.items():
		bamout = align_seq(args, sample, r, genome, outdir, 1)
		bamfiles1.append(bamout)
	# call the variants for the first time
	genome1 = call_snps(args, bamfiles1, genome, outdir, 1)
	
	# round2
	# prep sequence
        prepare_seq(args, genome1)
        # do the alignments for the first time
        bamfiles2 = []
        for sample, r in reads.items():
                bamout = align_seq(args, sample, r, genome1, outdir, 2)
                bamfiles2.append(bamout)
        # call the variants for the first time
        genome2 = call_snps(args, bamfiles2, genome1, outdir, 2)
	

        # round3
        # prep sequence
        prepare_seq(args, genome2)
        # do the alignments for the first time
        bamfiles3 = []
        for sample, r in reads.items():
                bamout = align_seq(args, sample, r, genome2, outdir, 3)
                bamfiles3.append(bamout)
        # call the variants for the first time
	genome3 = call_snps(args, bamfiles3, genome2, outdir, 3)

if __name__ == "__main__":
	main()
