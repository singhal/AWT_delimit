import re
import os
import pandas as pd
import glob
import subprocess

name = 'stacey_'
nodes = 1
cpu = 4
mem = 16
hours = 240

dir = '/scratch/drabosky_flux/sosi/AWT_delimit/stacey'

n = 0
dirs = os.listdir(dir)
for dir3 in dirs:
	newdir = os.path.join(dir, dir3)

	o = open("%s%s.sh" % (name, n), 'w')
	o.write("#PBS -N %s%s\n" % (name, n))
	o.write("#PBS -M sosi@umich.edu\n")
	o.write("#PBS -A drabosky_flux\n")
	o.write("#PBS -l qos=flux\n")
	o.write("#PBS -q flux\n")
	o.write("#PBS -l nodes=%s:ppn=%s,mem=%sgb\n" % (nodes, cpu, mem))
	o.write("#PBS -l walltime=%s:00:00\n" % hours)
	# o.write("#PBS -l walltime=06:14:21:58\n")
	o.write("#PBS -j oe\n")
	o.write("#PBS -V\n")
	o.write("cd %s\n" % newdir)

	files = glob.glob(newdir + '/*')
	o.write("~/bin/beast/bin/beast -threads 4 %s" % files[0])
	o.close()

	# if n > 1:
	#	subprocess.call("qsub %s%s.sh" % (name, n), shell=True)
	n += 1
