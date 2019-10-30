#!/usr/bin/python
import sys
import statistics
import numpy as np
import scipy.integrate as integrate
from scipy.integrate import quad

inFile = open(sys.argv[1],'r')
gene = inFile.readline()
num_windows = 0

def normalProbabilityDensity(x):
	constant = 1.0 / np.sqrt(2*np.pi)
	return(constant * np.exp((-x**2) / 2.0) )

f1_list = []
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
		if (all_f ==0):
			continue
		f1_list.append(f1_prob/all_f)
	gene = inFile.readline()
inFile.close()

# now convert to Z-score and sigma
global_mean_f1 = statistics.mean(f1_list)
global_std_f1 = statistics.stdev(f1_list)
print("F1 mean with statistics is " + str(global_mean_f1))
print("F1 pstdev with statistics is " + str(global_std_f1))


outFile = open("output/all_jan2014_25_to_32_starts_simulated_shifts_approach_4_theoretic_distributions", "w")
	
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

	# now go throught the simulated frameshifts, and compare the distributions of probabilities for them vs the un-augmented regions
	inFile = open("output/all_jan2014_25_to_32_starts_with_simulated_shifts_at_" + str(level), "r")

	gene = inFile.readline().strip()
	probs = []
	unsim_probs = []
	while gene:
		seq = inFile.readline()
		starts_string = inFile.readline()
		current_gene = gene
		gene = inFile.readline().strip()
		ss = [int(x) for x in starts_string.strip()[:-1].split(",")]
		if (current_gene in shift_info_dict[level]):
			# get the f1 fequency in start-stop region	
			f1_list = []
			f1_sum = 0
			total_sum = 0
			start = shift_info_dict[level][current_gene]["shift_start"]
			stop = shift_info_dict[level][current_gene]["shift_end"]
			for i in range(start, stop - 3, 3):
				f1_prob = int(ss[i])
				f2_prob = int(ss[i+1])
				f3_prob = int(ss[i+2])
				all_f = f1_prob + f2_prob + f3_prob
				if (all_f ==0):
					continue
				f1_sum += f1_prob
				total_sum += all_f
				f1_list.append(f1_prob/all_f)
			if (len(f1_list) == 0):
				#print("hit on an all-zero read region, no information")
				continue
			mean_f1 = statistics.mean(f1_list)
			#avg = f1_sum/total_sum # maybe later use avg as Steve suggested
			#print(current_gene + "\t" + str(mean_f1) + "\t" + str(avg) + "\n")
				
			# now get the z-score for the local mean	
			z_score = (mean_f1 - global_mean_f1)/global_std_f1
			prob = quad(normalProbabilityDensity, np.NINF, z_score) # the second argument is the error
			#print(current_gene + "\t" + str(mean_f1) + "\t" + str(z_score) + "\t" + str(prob[0]))
			probs.append(round(prob[0],3))

			if (level == 0.05):
				# now repeat everything for the unsimulated region of the level is 0.05
				ss = unsimulated_starts_dict[current_gene]
				# get the f1 fequency in start-stop region	
				f1_list = []
				f1_sum = 0
				total_sum = 0
				start = shift_info_dict[level][current_gene]["shift_start"]
				stop = shift_info_dict[level][current_gene]["shift_end"]
				for i in range(start, stop - 3, 3):
					f1_prob = int(ss[i])
					f2_prob = int(ss[i+1])
					f3_prob = int(ss[i+2])
					all_f = f1_prob + f2_prob + f3_prob
					if (all_f ==0):
						continue
					f1_sum += f1_prob
					total_sum += all_f
					f1_list.append(f1_prob/all_f)
				if (len(f1_list) == 0):
					#print("hit on an all-zero read region, no information")
					continue
				mean_f1 = statistics.mean(f1_list)
				#avg = f1_sum/total_sum # maybe later use avg as Steve suggested
				#print(current_gene + "\t" + str(mean_f1) + "\t" + str(avg) + "\n")
					
				# now get the z-score for the local mean	
				z_score = (mean_f1 - global_mean_f1)/global_std_f1
				prob = quad(normalProbabilityDensity, np.NINF, z_score) # the second argument is the error
				#print(current_gene + "\t" + str(mean_f1) + "\t" + str(z_score) + "\t" + str(prob[0]))
				unsim_probs.append(round(prob[0],3))
	if (level == 0.05):
		outFile.write(str(0.0) + "\n")
		outFile.write(",".join([str(x) for x in unsim_probs]) + "\n")
	outFile.write(str(level) + "\n")
	outFile.write(",".join([str(x) for x in probs]) + "\n")
outFile.close()
