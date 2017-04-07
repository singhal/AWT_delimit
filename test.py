import pandas as pd
import subprocess
import os
import glob

nodes = 1
cpu = 2
mem = 8
name = 'beast'
hours = 300

ix = 0
for num in [1,2,3]:
	for sp in ['robertsi', 'coggeri', 'carlia']:
		dirname = '/scratch/drabosky_flux/sosi/AWT_delimit/%s_phylogeny/starbeast/run%s/' % (sp, num)
		
        	sh_out = '%s%s.sh' % (name, ix)
        	o = open(sh_out, 'w')

        	o.write("#PBS -N %s%s\n" % (name, ix))
        	o.write("#PBS -M sosi@umich.edu\n")
        	o.write("#PBS -A drabosky_flux\n")
        	o.write("#PBS -l qos=flux\n")
        	o.write("#PBS -q flux\n")
        	o.write("#PBS -l nodes=%s:ppn=%s,mem=%sgb\n" % (nodes, cpu, mem))
        	o.write("#PBS -l walltime=%s:00:00\n" % hours)
        	o.write("#PBS -j oe\n")
        	o.write("#PBS -V\n")
        	o.write("\n")

		# o.write("module load python-anaconda2/latest\n")
		#  o.write("module load raxml/8.2.8\n")

		o.write("cd %s\n" % dirname)
		o.write("java -Xmx8g -jar /home/sosi/bin/beast/lib/beast.jar %s%s%s.xml\n" % (dirname, sp, num))
		o.close()
		ix += 1
