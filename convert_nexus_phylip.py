import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile", required=True)
args = parser.parse_args()

seq = {}
f = open(args.infile, 'r')
read = False
for l in f:
	if re.search('Matrix', l):
		read = True
	elif re.search('^;$', l):
		read = False
	else:
		if read:
			d = re.split('\s+', l.rstrip())
			seq[d[0]] = d[1]
f.close()

outfile = re.sub('nex', 'phy', args.infile)
o = open(outfile, 'w')
o.write('   %s %s\n' % (len(seq), len(list(seq.values())[0])))
for id, s in seq.items():
	o.write('%s\t%s\n' % (id, s))
o.close()
