#!/usr/bin/python
import sys
import statistics
import numpy as np
import scipy.integrate as integrate
from scipy.special import binom
import math

inFile = open(sys.argv[1],'r')
gene = inFile.readline()
num_windows = 0

p = 0.283
q = 0.717

#p = 0.05
#q = 0.95

# returns the exponent in scientific location the probablity - as the numbers are too small
def get_region_probability(starts, start, stop, gene="bla"):
	# correction for overly sharp spikes
	f23_sum = 0
	total_sum = 0
	for i in range(start, stop - 3, 3):
		f1_prob = int(starts[i])
		f2_prob = int(starts[i+1])
		f3_prob = int(starts[i+2])
		all_f = f1_prob + f2_prob + f3_prob
		total_sum += all_f
	for i in range(start, stop - 3, 3):
		if (int(starts[i+1]) >= total_sum/10):
			starts[i+1] = 0                       # zero out the super tall spikes
		if (int(starts[i+2]) >= total_sum/10):
			starts[i+2] = 0                       # zero out the super tall spikes

	# get the f23 counts in start-stop region	
	f23_sum = 0
	total_sum = 0
	for i in range(start, stop - 3, 3):
		f1_prob = int(starts[i])
		f2_prob = int(starts[i+1])
		f3_prob = int(starts[i+2])
		all_f = f1_prob + f2_prob + f3_prob
		if (all_f ==0):
			continue
		f23_sum += f2_prob + f3_prob
		total_sum += all_f
	if (total_sum == 0 or f23_sum == 0):
		#print("hit on an all-zero read region, no information")
		return "UNDEF",total_sum,0
	#print(gene)	
	#print(str(total_sum) + "\t" + str(f23_sum))
	#  try 1 : binomial formula - overflow
	#for i in range(f23_sum, total_sum+1):
	#	res = binom(total_sum, i) * pow(p, i) * pow(q, total_sum-i)
	#	sum_prob += res 

	# try 2: logs in sum - too slow
	#total_prob = 0	
	#for i in range(f23_sum, total_sum):
	#	sum_prob = 0
	#	for j in range(1, total_sum + 1):
	#		sum_prob += math.log10(j)
	#	for j in range(1, total_sum-i+1):
	#		sum_prob -= math.log10(j)
	#	for j in range(1, i+1):
	#		sum_prob -= math.log10(j)
	#	sum_prob += i * math.log10(p)
	#	sum_prob += (total_sum - i) * math.log10(q)
	#	total_prob += math.exp(sum_prob)
	#print(total_prob)

	# try 3: stirling forumula, fine but neeed to get rid of the log 
	# take math.log10s of binomial equation
	#total_prob = 0	
	#for i in range(f23_sum, total_sum):
	#	sum_prob = 0
	#	sum_prob += total_sum * math.log10(total_sum) - total_sum + 1
	#	sum_prob -= (total_sum - i) * math.log10(total_sum - i) - (total_sum - i) + 1
	#	sum_prob -= i * math.log10(i) - i + 1
	#	sum_prob += i * math.log10(p)
	#	sum_prob += (total_sum - i) * math.log10(q)
	#	total_prob += sum_prob


	# try 4 - now try to get just the upper limit on the probability - the square
	total_prob = 0	
	total_prob += total_sum * math.log10(total_sum) - total_sum + 1
	total_prob -= (total_sum - f23_sum) * math.log10(total_sum - f23_sum) - (total_sum - f23_sum) + 1
	total_prob -= f23_sum * math.log10(f23_sum) - f23_sum + 1
	total_prob += f23_sum * math.log10(p)
	total_prob += (total_sum - f23_sum) * math.log10(q)
	# now multiply by the width of the box - equivalent to adding the exponents
	total_prob += math.log10(total_sum - f23_sum)
	total_p = f23_sum/total_sum
	#if (gene == "YGR192C"):
	#	print(str(start) + "\t" + str(stop))
	#	print(str(total_sum) + "\t" + str(f23_sum))
	#	print(total_prob)
	return min(0, total_prob),total_sum,total_p

outFile = open("output/all_jan2014_25_to_32_starts_simulated_shifts_approach_5_theoretic_distributions_with_before", "w") #_higherf1
outFile2 = open("output/all_jan2014_25_to_32_starts_simulated_shifts_approach_5_averages_of_distributions_with_before", "w") #_higherf1
	
# save the unsimulated information
unsimulated_starts_dict = dict()
inFile = open("output/all_jan2014_25_to_32_starts", "r")
gene = inFile.readline().strip()
while gene:
	seq = inFile.readline()
	starts_string = inFile.readline()
	current_gene = gene
	gene = inFile.readline().strip()
	ss = [int(x) for x in starts_string.strip()[:-1].split(",")]
	unsimulated_starts_dict[current_gene] = ss

# get the shifts information
shift_info_dict = dict()
levels = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
#levels = [0.05]
for level in levels:
	infoFile = open("output/all_jan2014_25_to_32_starts_simulated_shifts_info_at_" + str(level), "r")
	shift_info_dict[level] = dict()
	for line in infoFile:
		cs = line.strip().split(",")
		gene = cs[0]
		shift_start = int(cs[1])
		shift_end = int(cs[2])
		direction = int(cs[3])
		percent = float(cs[4])
		#if (percent == '0.9'):
		shift_info_dict[level][gene] = dict()
		shift_info_dict[level][gene]["shift_start"] = shift_start
		shift_info_dict[level][gene]["shift_end"] = shift_end
		shift_info_dict[level][gene]["direction"] = direction
		shift_info_dict[level][gene]["percent"] = percent

	# now go throught the simulated frameshifts, and do the binomial probabilities
	inFile = open("output/all_jan2014_25_to_32_starts_with_simulated_shifts_at_" + str(level), "r")

	gene = inFile.readline().strip()
	probs = []
	unsim_probs = []
	xs = []
	ys = []
	colors = []
	unsim_xs = []
	unsim_ys = []
	while gene:
		seq = inFile.readline()
		starts_string = inFile.readline()
		current_gene = gene
		gene = inFile.readline().strip()
		ss = [int(x) for x in starts_string.strip()[:-1].split(",")]
		if (current_gene in shift_info_dict[level]):
			# get the f23 counts in start-stop region	
			start = shift_info_dict[level][current_gene]["shift_start"]
			stop = shift_info_dict[level][current_gene]["shift_end"]

			frameshift_prob,total_sum,total_p = get_region_probability(ss, start, stop)
			if (frameshift_prob == "UNDEF"):
				continue
		
			before_val = max(0, start - (stop - start)) - max(0, start - (stop - start))%3	
			before_frameshift_prob,before_sum,before_p = get_region_probability(ss, before_val, start)
			if (before_frameshift_prob == "UNDEF"):
				continue

			total_prob = min(0, (frameshift_prob - before_frameshift_prob))
			# don't consider situtations where 23 counts are not greater than average, and where
			if ((before_p > total_p) or (total_p < p)):
				total_prob = 0
			
			probs.append(total_prob)
			xs.append(total_sum)
			ys.append(total_prob)
			if ((stop-start) < 60):
				colors.append("yellow")
			elif ((stop-start) < 80):
				colors.append("green")
			elif ((stop-start) < 100):
				colors.append("blue")
			else:
				colors.append("red")
		

			if (level == 0.05):
				# now repeat everything for the unsimulated region of the level is 0.05
				ss = unsimulated_starts_dict[current_gene]

				start = shift_info_dict[level][current_gene]["shift_start"]
				stop = shift_info_dict[level][current_gene]["shift_end"]

				frameshift_prob,total_sum,total_p = get_region_probability(ss, start, stop, current_gene)
				if (frameshift_prob == "UNDEF"):
					continue
			
				before_val = max(0, start - (stop - start)) - max(0, start - (stop - start))%3	
				before_frameshift_prob,before_sum,before_p = get_region_probability(ss, before_val, start, current_gene)
				if (before_frameshift_prob == "UNDEF"):
					continue

				total_prob = min((frameshift_prob -  before_frameshift_prob), 0)
				# don't consider situtations where 23 counts are not greater than average, and where
				if ((before_p > total_p) or (total_p < p)):
					total_prob = 0

				if (total_prob <= -30):
					print(current_gene)
					print(total_prob)
				unsim_probs.append(total_prob)
				unsim_xs.append(total_sum)
				unsim_ys.append(total_prob)
	if (level == 0.05):
		outFile.write(str(0.0) + "\n")
		outFile.write(",".join([str(x) for x in unsim_probs]) + "\n")
		outFile2.write(str(0.0) + "\n")
		outFile2.write(",".join([str(x) for x in unsim_xs]) + "\n")
		outFile2.write(",".join([str(x) for x in unsim_ys]) + "\n")
		outFile2.write(",".join(colors) + "\n")
	outFile.write(str(level) + "\n")
	outFile.write(",".join([str(x) for x in probs]) + "\n")
	
	# calculate the average_bins
	outFile2.write(str(level) + "\n")
	outFile2.write(",".join([str(x) for x in xs]) + "\n")
	outFile2.write(",".join([str(x) for x in ys]) + "\n")
	outFile2.write(",".join(colors) + "\n")
		
outFile.close()
