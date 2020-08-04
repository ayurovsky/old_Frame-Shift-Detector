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
import scipy.integrate as integrate
from scipy.special import binom
import math
sys.stdout.flush()

# for 4 good datasets, 28 and 29 reads only
p = 0.086
q = 0.914
p_f2 = 0.045
p_f3 = 0.041
# so, the range

# for all_jan, 28 reads only
#p = 0.09
#q = 0.91

# for Wu2019, and Guydhosh2014, only 28 and 29 reads
#p = 0.08
#q = 0.92

# for Wu2019, and Guydhosh2014, only 28 and 29 reads
#p = 0.08
#q = 0.92

# for all seven datasets combined
#p = 0.241 
#q = 0.759

# for all but wu and jan
#p = 0.239 
#q = 0.761


# for Young2015
#p = 0.212
#q = 0.788

# for Williams2014 
#p = 0.238 
#q = 0.762

# for Hanson2018 , and Geraschenko2014
#p = 0.298
#q = 0.702

# for allJan2016
#p = 0.283
#q = 0.717

p_value_threshold = -5# -7 #-4

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
	f2_sum = 0
	f3_sum = 0
	for i in range(start, stop - 3, 3):
		f1_prob = int(starts[i])
		f2_prob = int(starts[i+1])
		f3_prob = int(starts[i+2])
		all_f = f1_prob + f2_prob + f3_prob
		if (all_f ==0):
			continue
		f23_sum += f2_prob + f3_prob
		f2_sum += f2_prob
		f3_sum += f3_prob
		total_sum += all_f
	if (total_sum == 0 or f23_sum == 0 or f23_sum == total_sum):
		#print("hit on an all-zero read region, no information")
		return "UNDEF",total_sum,0,0,0
	#  try 1 : binomial formula - overflow
	# try 2: logs in sum - too slow
	# try 3: stirling forumula, fine but neeed to get rid of the log 
	# try 4 - now try to get just the upper limit on the probability - the square
	total_prob = 0	
#	print(total_sum)
#	print(f23_sum)
#	print("\n")
	total_prob += total_sum * math.log10(total_sum) - total_sum + 1
	total_prob -= (total_sum - f23_sum) * math.log10(total_sum - f23_sum) - (total_sum - f23_sum) + 1
	total_prob -= f23_sum * math.log10(f23_sum) - f23_sum + 1
	total_prob += f23_sum * math.log10(p)
	total_prob += (total_sum - f23_sum) * math.log10(q)
	# now multiply by the width of the box - equivalent to adding the exponents
	total_prob += math.log10(total_sum - f23_sum)
	total_p = f23_sum/total_sum
	return min(0, total_prob),total_sum,total_p,f2_sum,f3_sum

def find_stop_codon(sequence, start, end):
	stop = -1
	for i in range(start,end,3):
		if (sequence[i:i+3] in ["TGA", "TAA", "TAG"]):
			stop = i
			break
	return stop

inFile = open(sys.argv[1],'r')
outFile = open(sys.argv[1] + "_with_found_frameshifts_approach_5.5", 'w')

gene = inFile.readline().strip()
counter = 0
while gene:
	seq = inFile.readline()
	starts_string = inFile.readline()
	current_gene = gene
	gene = inFile.readline().strip()

	ss = starts_string.strip()[:-1].split(",")
	best_prob = 0 
	best_info = "" 

	for i in range(0, len(ss) - 3, 3): # i is in frame

		start = i

		# check 2nd frame for the +1 shift
		stop = find_stop_codon(seq, start+1, len(seq))
		plus_1_length = stop - start 
		frameshift_prob,total_sum_plus_1,total_p,f2_sum_for_plus_1,f3_sum_for_plus_1 = get_region_probability(ss, start, stop)
		if (total_sum_plus_1 == 0):
			percent_shifted_f2 = 0	
		else:
			percent_shifted_f2 = min(1.0, (f2_sum_for_plus_1/total_sum_plus_1 - p_f2)/(q - p_f2))
	
		before_val = max(0, start - (stop - start)) - max(0, start - (stop - start))%3	
		before_frameshift_prob,before_sum_plus_1,before_p,f2_sum,f3_sum = get_region_probability(ss, before_val, start)

		total_prob = 0 
		# try to assign value to the exponent of probability
		if ((before_frameshift_prob != "UNDEF") and (frameshift_prob != "UNDEF") and ((stop - start) > 30)):
			total_prob = min(0, (frameshift_prob - before_frameshift_prob))
			# don't consider situtations where 23 counts are not greater than average, and where before counts are greater than after counts
			# also where there are more reads in a the wrong frame
			if ((before_p > total_p) or (total_p < p) or (f2_sum_for_plus_1 < f3_sum_for_plus_1)):
				total_prob = 0
		plus_1_prob = total_prob
		before_plus_1_prob = before_frameshift_prob

		# check the 3rd frame for the -1 shift	
		stop = find_stop_codon(seq, start-1, len(seq))
		minus_1_length = stop - start 
		frameshift_prob,total_sum_minus_1,total_p,f2_sum_minus_1,f3_sum_minus_1 = get_region_probability(ss, start, stop)
		if (total_sum_minus_1 == 0):
			percent_shifted_f3 = 0	
		else:
			percent_shifted_f3 = min(1.0, (f3_sum_minus_1/total_sum_minus_1 - p_f3)/(q - p_f3))
	
		before_val = max(0, start - (stop - start)) - max(0, start - (stop - start))%3	
		before_frameshift_prob,before_sum_minus_1,before_p,f2_sum,f3_sum = get_region_probability(ss, before_val, start)

		total_prob = 0 
		# try to assign value to the exponent of probability
		if ((before_frameshift_prob != "UNDEF") and (frameshift_prob != "UNDEF") and ((stop - start) > 30)):
			total_prob = min(0, (frameshift_prob - before_frameshift_prob))
			# don't consider situtations where 23 counts are not greater than average, and where before counts are greater than after counts
			# also where there are more reads in a the wrong frame
			if ((before_p > total_p) or (total_p < p) or (f3_sum_minus_1 < f2_sum_minus_1)):
				total_prob = 0
		minus_1_prob = total_prob
		before_minus_1_prob = before_frameshift_prob

		# now see if we are doing better than the previous attempts
		real_length = plus_1_length
		direction = "plus"
		percent_shifted = percent_shifted_f2
		total_sum = total_sum_plus_1
		before_prob = before_plus_1_prob
		if (minus_1_prob < plus_1_prob):
			real_length = minus_1_length
			direction = "minus"
			percent_shifted = percent_shifted_f3
			total_sum = total_sum_minus_1
			before_prob = before_minus_1_prob
			
		final_prob = min(plus_1_prob, minus_1_prob)	
		ustart = "no"
		if (start < 30):
			ustart = "yes"

		if (final_prob <= p_value_threshold): #-4):
			# fake the info, don't care now, using just for plots
			outString = ""
			outString += "Gene: " + current_gene + "\tp-value: " + str(final_prob) + "\tgene-length: " + str(len(seq)) + "\tshift-start-position: " + str(start) + "\tshift-length " + str(real_length) + "\tdirection " + direction + "\tbefore_gene_start " + ustart + "\tpercent_shifted " + str(percent_shifted) + "\tnum_reads " + str(total_sum) + "\tbefore_region_p-value " + str(before_prob) + "\n"
			outString += "seq_info" + "\n"
			outString += "f1_starts" + "\n"
			outString += "f_alt_starts" + "\n"
			outString += "percents" + "\n"
			outString += "\n"

			if (final_prob < best_prob):
				best_prob = final_prob
				best_info = outString
	if (best_prob <= p_value_threshold): #-4):
		outFile.write(best_info)
	counter += 1
	#if (counter %100 == 0):
	#	print(counter)
inFile.close()
outFile.close()
