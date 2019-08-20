import subprocess
import os
from subprocess import Popen, PIPE
import sys

number_genes = 4612
genes_per_job = 10

#max_jobs = 500
processes = []

for m in range(0, int(number_genes/genes_per_job)+1):
	start = m * genes_per_job 
	end = start + genes_per_job 
	#args = ['python', 'find_frameshifts.py', 'bla_28', str(start), str(end)]
	args = ['python', 'find_frameshifts.py', 'bla_25_to_32', str(start), str(end)]
	p = subprocess.Popen(args, stdout=PIPE, stderr=PIPE)
	processes.append(p)
	#print(args)

for p in processes:
	stdout, stderr = p.communicate()
print("done")
