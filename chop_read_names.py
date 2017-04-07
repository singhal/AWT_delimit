import re
import os
import gzip
import subprocess

sample = 'SS136'
dir = '/scratch/drabosky_flux/sosi/AWT_delimit/trim_reads/'

def clean_file(sample, end1, end2):
	file1 = os.path.join(dir, '%s_R%s.final.fq.gz' % (sample, end1))
	out1 = os.path.join(dir, '%s_tmp%s.fq.gz' % (sample, end1))

	i = gzip.open(file1, 'r')
	o = gzip.open(out1, 'w')

	while True:
		line1 = i.readline()

		if not line1: break

		line2 = i.readline().rstrip()
                line3 = i.readline().rstrip()
                line4 = i.readline().rstrip()

		line1 = re.sub('_%s' % end2, '', line1.rstrip())

		o.write('%s\n%s\n%s\n%s\n' % (line1, line2, line3, line4))

	i.close()
	o.close()

	subprocess.call("mv %s %s" % (out1, file1), shell=True)

clean_file(sample, 1, 'F')
clean_file(sample, 2, 'R')

		
		
