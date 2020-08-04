#!/usr/bin/python
import sys
import statistics
import scipy.integrate as integrate
from scipy.integrate import quad

inFile = open(sys.argv[1],'r')
gene = inFile.readline()


f1_list = []
total_f1 = 0
total_f2 = 0
total = 0
# get the f1 now
while gene:
	seq = inFile.readline()
	starts_string = inFile.readline()
	ss = starts_string.split(",")
	
	for i in range(0, len(seq) - 3, 3): # i is in frame
		f1_prob = int(ss[i])
		f2_prob = int(ss[i+1])
		f3_prob = int(ss[i+2])
		all_f = f1_prob + f2_prob + f3_prob
		total_f1 += f1_prob
		total_f2 += f2_prob
		total += all_f
		if (all_f ==0):
			continue
		f1_list.append(f1_prob/all_f)
	gene = inFile.readline()
inFile.close()

global_mean_f1 = statistics.mean(f1_list)
global_std_f1 = statistics.stdev(f1_list)
print("F1 mean with statistics is " + str(global_mean_f1))
print("F1 stdev with statistics is " + str(global_std_f1))
print("Total f1 is " + str(total_f1/total))
print("Total f2 is " + str(total_f2/total))
print("Total is: " + str(total))
