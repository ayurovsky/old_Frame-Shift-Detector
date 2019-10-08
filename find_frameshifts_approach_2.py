#!/usr/bin/python
import sys
import random
from scipy import stats 
import numpy as np
from statistics import mean
sys.stdout.flush()

def get_drop(loc, starts): 
	# running average of 11 at loc 
	f1 = 0
	f2 = 0
	f3 = 0
	for j in range(max(0,loc-15),min(len(starts)-3,loc+15),3):
		f1 += int(starts[j])
		f2 += int(starts[j+1])
		f3 += int(starts[j+2])
	all_f = f1 + f2 + f3
	loc_percent = 0.0
	if (all_f !=0):
		loc_percent = f1/all_f
	# running average of 11 before loc 
	f1 = 0
	f2 = 0
	f3 = 0
	for j in range(max(0,loc-1-15),min(len(starts)-3,loc-1+15),3):
		f1 += int(starts[j])
		f2 += int(starts[j+1])
		f3 += int(starts[j+2])
	all_f = f1 + f2 + f3
	before_percent = 0.0
	if (all_f !=0):
		before_percent = f1/all_f

	drop = round(before_percent - loc_percent, 3)
	return drop

def find_stop_codon(sequence, start, end):
	stop = -1
	for i in range(start,end,3):
		if (sequence[i:i+3] in ["TGA", "TAA", "TAG"]):
			stop = i
			break
	return stop

def custom_get_percent_in_range(starts, start, end, frame="f1"):
	f1 = 0
	f2 = 0
	f3 = 0
	for j in range(start, end, 3): # i is in frame
		f1 += int(starts[j])
		f2 += int(starts[j+1])
		f3 += int(starts[j+2])
	all_f = f1 + f2 + f3
	if (all_f !=0):
		if (frame == "f1"):
			return f1/all_f
		if (frame == "f2"):
			return f2/all_f
		if (frame == "f3"):
			return f3/all_f
	else:
		return 0.0

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

inFile = open(sys.argv[1],'r')
gene = inFile.readline().strip()
outFile = open(sys.argv[1] + "_with_found_frameshifts_approach_2", 'w')
badFile = open(sys.argv[1] + "_without_threshold_frameshifts_approach_2", 'w')

counter = 0
while gene:
	seq = inFile.readline()
	starts_string = inFile.readline()
	current_gene = gene
	gene = inFile.readline().strip()
	#if (current_gene != "YBL041W"):
	#	continue

	ss = starts_string.strip()[:-1].split(",")

	found_frameshifts_counter = 0

	best_p_value = 1
	best_frameshift_info = ""

	for i in range(0, len(ss) - 3, 3): # i is in frame

		# get drop value in this location 
		drop_value = get_drop(i, ss)
		#print(str(i) + ": " + str(drop_value))

		if (drop_value > 0): # there is a decrease

			drop = i	
			before_percents = np.array(custom_get_range_of_percents(ss, 0, drop, "f1"))

			# decide if this is a + 1 or a -1 shift

			# check 2nd frame for the +1 shift
			plus_1_stop = find_stop_codon(seq, drop+1, len(seq))
			after_percents = np.array(custom_get_range_of_percents(ss, drop, plus_1_stop, "f1"))
			t_test = stats.ttest_ind(before_percents, after_percents, equal_var=False)
			p_value_plus_1 = t_test[1]


			# check 3rd frame for -1 shift
			minus_1_stop = find_stop_codon(seq, drop-1, len(seq))
			after_percents = np.array(custom_get_range_of_percents(ss, drop, minus_1_stop, "f1"))
			t_test = stats.ttest_ind(before_percents, after_percents, equal_var=False)
			p_value_minus_1 = t_test[1]

			p_value = 1
			if ((p_value_plus_1 < p_value_minus_1) and ((plus_1_stop - drop) > 30)):
				p_value = p_value_plus_1
				real_new_frame = "f2"
				real_length = plus_1_stop - drop
			elif ((p_value_minus_1 < p_value_plus_1) and ((minus_1_stop - drop) > 30)):
				p_value = p_value_minus_1
				real_new_frame = "f3"
				real_length = minus_1_stop - drop
			
				
			if (p_value <= 0.0001):
				found_frameshifts_counter += 1			

				outString = ""
				outString += "Gene: " + current_gene + "\tp-value: " + str(p_value) + "\tgene-length: " + str(len(seq)) + "\tshift-start-position: " + str(drop) + "\tshift-length " + str(real_length) + "\n"
				f1_starts = "" 
				f_alt_starts = ""
				seq_info = ""
				for i in range(drop-30,drop+30,3):
					if (i == drop):
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
				
				if (p_value < best_p_value):
					best_p_value = p_value
					best_frameshift_info = outString
	#print(current_gene + " has " + str(found_frameshifts_counter) + " frameshifts above threshold\n")
	if (best_p_value < 1):
		outFile.write(best_frameshift_info)
	else:
		outString = "Gene: " + current_gene + "\n"
		percent_f1_all = ",".join(str(round(x,3)) for x in custom_get_range_of_percents(ss, 0, len(ss)-2, "f1"))
		outString += percent_f1_all
		outString += "\n\n"
		badFile.write(outString)
inFile.close()
badFile.close()
outFile.close()
