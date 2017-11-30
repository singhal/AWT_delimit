import re
import os
import pandas as pd
import glob
import subprocess

name = 'count_'
nodes = 1
cpu = 1
mem = 4
hours = 240

files = glob.glob('/scratch/drabosky_flux/sosi/AWT_delimit/variants/*gz')
files = [file for file in files if re.search('qual_filtered.cov_filtered.vcf.gz', file)] 

n = 0
for file in files:
        # newdir = os.path.join(dir, dir3)

        o = open("%s%s.sh" % (name, n), 'w')
        o.write("#PBS -N %s%s\n" % (name, n))
        o.write("#PBS -M sosi@umich.edu\n")
        o.write("#PBS -A drabosky_flux\n")
        o.write("#PBS -l qos=flux\n")
        o.write("#PBS -q flux\n")
        o.write("#PBS -l nodes=%s:ppn=%s,mem=%sgb\n" % (nodes, cpu, mem))
        o.write("#PBS -l walltime=%s:00:00\n" % hours)
        o.write("#PBS -j oe\n")
        o.write("#PBS -V\n")
        o.write("python ~/AWT_delimit/count_depth.py --file %s\n" % file)
	n += 1
