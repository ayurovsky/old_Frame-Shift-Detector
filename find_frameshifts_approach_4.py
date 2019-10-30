#!/usr/bin/python
import sys
import random
from scipy import stats 
import numpy as np
from statistics import mean
import statistics
import numpy as np
import scipy.integrate as integrate
from scipy.integrate import quad
sys.stdout.flush()

global_mean_f1 = 0
global_std_f1 = 0

def custom_get_range_of_percents(starts, start, end, frame="f1"):
	percents = []
	for j in range(start, end, 3): # i is in frame
		f1 = int(starts[j])
		f2 = int(starts[j+1])
		f3 = int(starts[j+2])
		all_f = f1 + f2 + f3
		if (all_f !=0):
			if (frame == "f1"):
				percents.append(f1/all_f)
			if (frame == "f2"):
				percents.append(f2/all_f)
			if (frame == "f3"):
				percents.append(f3/all_f)
		else:
			percents.append(0.0)
	return percents

def normalProbabilityDensity(x):
	constant = 1.0 / np.sqrt(2*np.pi)
	return(constant * np.exp((-x**2) / 2.0) )

def find_stop_codon(sequence, start, end):
	stop = -1
	for i in range(start,end,3):
		if (sequence[i:i+3] in ["TGA", "TAA", "TAG"]):
			stop = i
			break
	return stop

def get_probability(starts, start, stop):
	f1_list = []
	f1_sum = 0
	for i in range(start, stop - 3, 3):
		f1_prob = int(starts[i])
		f2_prob = int(starts[i+1])
		f3_prob = int(starts[i+2])
		all_f = f1_prob + f2_prob + f3_prob
		if (all_f ==0):
			continue
		f1_sum += f1_prob
		f1_list.append(f1_prob/all_f)
	if (len(f1_list) == 0):
		return(0)	
	mean_f1 = statistics.mean(f1_list)
	# now get the z-score for the local mean	
	z_score = (mean_f1 - global_mean_f1)/global_std_f1
	prob = quad(normalProbabilityDensity, np.NINF, z_score) # the second argument is the error
	return(prob[0])

#####################################################################

# calculate the global mean and std f1
inFile = open(sys.argv[1],'r')

f1_list = []
gene = inFile.readline()
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

#####################################################################

genes_dict = dict()
infoFile = open("output/all_jan2014_25_to_32_starts_simulated_shifts_info_at_0.5", "r")
for line in infoFile:
	cs = line.strip().split(",")
	gene = cs[0]
	genes_dict[gene] = "yes"

inFile = open(sys.argv[2],'r')
outFile = open(sys.argv[2] + "_with_found_frameshifts_approach_4", 'w')

counter = 0
gene = inFile.readline().strip()
while gene:
	seq = inFile.readline()
	starts_string = inFile.readline()
	current_gene = gene
	gene = inFile.readline().strip()

	ss = starts_string.strip()[:-1].split(",")

	found_frameshifts_counter = 0
	best_diff = 0 
	best_prob = 1
	best_pvalue = 1

	#if (current_gene not in genes_dict):
	#continue
	#if (current_gene != "YAL012W"):
	#	continue

	for i in range(0, len(ss) - 3, 3): # i is in frame

		# decide if this is a + 1 or a -1 shift

		# check 2nd frame for the +1 shift
		plus_1_stop = find_stop_codon(seq, i+1, len(seq))
		after_probability_plus_1 = get_probability(ss, i, plus_1_stop)
		plus_1_length = plus_1_stop - i
		before_val = max(0, i - plus_1_length) - max(0, i - plus_1_length)%3	
		before_probability_plus_1 = get_probability(ss, before_val, i) 

		# check 3rd frame for -1 shift
		minus_1_stop = find_stop_codon(seq, i-1, len(seq))
		after_probability_minus_1 = get_probability(ss, i, minus_1_stop)
		minus_1_length = minus_1_stop - i
		before_val = max(0, i - minus_1_length) - max(0, i - minus_1_length)%3	
		before_probability_minus_1 = get_probability(ss, before_val, i) 

		prob = 1 
		real_length = 0
		p_value = 1	
	
		if (((before_probability_plus_1 - after_probability_plus_1) > (before_probability_minus_1 - after_probability_minus_1)) and ((plus_1_stop - i) > 30)):
			prob = after_probability_plus_1
			real_new_frame = "f2"
			real_length = plus_1_length
			before_probability = before_probability_plus_1
			after_percents = np.array(custom_get_range_of_percents(ss, i, plus_1_stop, "f1"))
			new_length = int( 3 * round( (plus_1_stop - i) / 3. ))
			before_percents = np.array(custom_get_range_of_percents(ss, max(0, i - new_length), i, "f1"))
			t_test = stats.ttest_ind(before_percents, after_percents, equal_var=False)
			p_value = t_test[1]
		elif (((before_probability_plus_1 - after_probability_plus_1) < (before_probability_minus_1 - after_probability_minus_1)) and ((minus_1_stop - i) > 30)):
			prob = after_probability_minus_1
			real_new_frame = "f3"
			real_length = minus_1_length
			before_probability = before_probability_minus_1
			after_percents = np.array(custom_get_range_of_percents(ss, i, minus_1_stop, "f1"))
			new_length = int( 3 * round( (minus_1_stop - i) / 3. ))
			before_percents = np.array(custom_get_range_of_percents(ss, max(0, i - new_length), i, "f1"))
			t_test = stats.ttest_ind(before_percents, after_percents, equal_var=False)
			p_value = t_test[1]
		
		# check that prob is within limits - not zero (that means no reads), and less than limit	
		if (prob > 0 and prob <= 0.4):
			found_frameshifts_counter += 1			

			outString = ""
			outString += "Gene: " + current_gene + "\tp-value: " + str(p_value) + "\tgene-length: " + str(len(seq)) + "\tshift-start-position: " + str(i) + "\tshift-length " + str(real_length) + "\t" + real_new_frame + "\tprob: " + str(prob) + "\tbefore_prob: " + str(before_probability) + "\n"
			#outString += "Gene: " + current_gene + "\t" + real_new_frame + "\tprob: " + str(prob) + "\tbefore_prob: " + str(before_probability) + "\tgene-length: " + str(len(seq)) + "\tshift-start-position: " + str(i) + "\tshift-length " + str(real_length) + "\n"
			f1_starts = "" 
			f_alt_starts = ""
			seq_info = ""
			for i in range(i-30,i+30,3):
				if (i == i):
					f1_starts += "| "
					f_alt_starts += "| "
					seq_info += "| "
				f1_starts += str(ss[i]) + "   "
				seq_info += seq[i:i+1] + " " + seq[i+1:i+2] + " " + seq[i+2:i+3] + " "
				if (real_new_frame == "f2"):
					f_alt_starts += " " + str(ss[i+1]) + "  " 
				else:
					f_alt_starts += "  " + str(ss[i-1]) + " "
			outString += seq_info + "\n"
			outString += f1_starts + "\n"
			outString += f_alt_starts + "\n"
			percent_f1_all = ",".join(str(round(x,3)) for x in custom_get_range_of_percents(ss, 0, len(ss)-2, "f1"))
			outString += percent_f1_all
			outString += "\n\n"
			

			if ((before_probability - prob) > best_diff):
				best_diff = before_probability - prob
				best_frameshift_info = outString
				best_p_value = p_value
	if (best_diff >= 0.1 and best_p_value <= 0.0001):
		outFile.write(best_frameshift_info)
	counter += 1
	if (counter %100 == 0):
		print(counter)
inFile.close()
outFile.close()
