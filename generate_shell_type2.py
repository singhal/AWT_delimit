import pandas as pd
import subprocess
import os
import glob

nodes = 1
cpu = 4
mem = 16
name = 'phylo'
hours = 100

'''
cmds = []
f = open('commands', 'r')
for l in f:
	cmds.append(l.rstrip())
f.close()

for ix, cmd in enumerate(cmds):
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

	o.write(cmd)
	o.close()
'''
ix = 1
for num in range(0,20):
	for sp in ['robertsi', 'coggeri', 'carlia']:
		dirname = '/scratch/drabosky_flux/sosi/AWT_delimit/original_approach/%s/run%s/' % (sp, num)
		file = glob.glob(dirname + '*nex')[0]
		
        	sh_out = '%s%s.sh' % (name, ix)
        	o = open(sh_out, 'w')

        	o.write("#PBS -N %s%s\n" % (name, ix))
        	o.write("#PBS -M sosi@umich.edu\n")
        	o.write("#PBS -A sosi_flux\n")
        	o.write("#PBS -l qos=flux\n")
        	o.write("#PBS -q flux\n")
        	o.write("#PBS -l nodes=%s:ppn=%s,mem=%sgb\n" % (nodes, cpu, mem))
        	o.write("#PBS -l walltime=%s:00:00\n" % hours)
        	o.write("#PBS -j oe\n")
        	o.write("#PBS -V\n")
        	o.write("\n")

		#  o.write("module load python-anaconda2/latest\n")
		#  o.write("module load raxml/8.2.8\n")

		# o.write("cd %s\n" % dirname)
		# o.write("java -Xmx8g -jar /home/sosi/bin/beast/lib/beast.jar %s%s%s.xml\n" % (dirname, sp, num))
		o.write("module load  beagle/201606\n")
		o.write("java -Xmx16G -jar $PHYLONET ~/bin/PhyloNet_3.6.1.jar %s" % (file))
		# o.close()
		ix += 1

		subprocess.call("qsub %s" % sh_out, shell=True)
