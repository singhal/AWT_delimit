import re
import gzip
import os
import argparse

# number of loci	number of high-coverage sites	average coverage

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", required=True, help='File with the sample data.')
args = parser.parse_args()

infile = args.file
f = gzip.open(infile, 'r')

d = {}
for l in f:
	if re.search('^#CHROM', l):
		inds = re.split('\t', l.rstrip())
		inds = inds[9:]
		
		for ind in inds:
			d[ind] = {'loci': {}, 'num_sites': 0, 'coverage': 0}
	elif not re.search('^#', l):
		x = re.split('\t', l.rstrip())
		data = re.split(':', x[8])
		dp = data.index('DP')

		for geno, ind in zip(x[9:], inds):
			data = re.split(':', geno)
			if not data[0] == './.':
				# print(data, dp)
				depth = int(data[dp])
	
				if depth >= 10:
					d[ind]['num_sites'] += 1
					d[ind]['coverage'] += depth
					d[ind]['loci'][x[0]] = 1
f.close()

for ind in d:
	cov = d[ind]['coverage'] / float( d[ind]['num_sites'] )
	print('%s,%s,%s,%s' % (ind, len(d[ind]['loci']), d[ind]['num_sites'], cov))
			

