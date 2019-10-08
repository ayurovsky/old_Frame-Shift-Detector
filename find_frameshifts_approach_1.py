#!/usr/bin/python
import sys
import random
from scipy import stats 
import numpy as np
from statistics import mean

def get_drop(sequence, starts, print_output=False): 
	print_string = ""
	values_string = []
	for i in range(0, len(starts) - 3, 3): # i is in frame
		# for running average of eleven 
		f1 = 0
		f2 = 0
		f3 = 0
		for j in range(max(0,i-15),min(len(starts)-3,i+15),3):
			f1 += int(starts[j])
			f2 += int(starts[j+1])
			f3 += int(starts[j+2])
		all_f = f1 + f2 + f3
		percent = 0.0
		if (all_f !=0):
			percent = f1/all_f
		print_string += str(round(percent,3)) + ","	
		values_string.append(percent)
	if (print_output):
		print("Percent running average of 11AA")
		print(print_string)
	drops_string = [0]
	max_i = -1
	max_value = -1
	for i in range(1, len(values_string)): # i is in frame
		drop = values_string[i-1] - values_string[i]
		drops_string.append(round(drop,3))
		if (drop > max_value and values_string[i] != 0.0):
			max_value = drop
			max_i = i
	max_i *= 3 # back into the sequence space
	return (max_i,max_value,values_string) 

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
#start_idx = int(sys.argv[2])
#end_idx = int(sys.argv[3])
outFile = open(sys.argv[1] + "_with_found_frameshifts_approach_1", 'w')

counter = 0
while gene:
	seq = inFile.readline()
	starts_string = inFile.readline()
	current_gene = gene
	gene = inFile.readline().strip()
	#if (not (current_gene == "YIL009C-A" or current_gene == "YOR239W")):
	#	continue
	#if (counter < start_idx or counter >= end_idx):
	#	counter += 1
	#	continue
	#else:
	#	counter += 1

	ss = starts_string.strip()[:-1].split(",")

	# get the max drop - % change in frame 1 reads	
	#print("_" + current_gene + "_")
	#print(seq)
	drop_info = get_drop(seq, ss)
	#if (current_gene == "YDL137W"):
	#	drop_info = get_drop(seq, ss, True)
	#	print("Max drop location for running 11 avearge of percent location is: " + str(drop_info[0]) + " " + str(drop_info[1]))

	# decide if this is a + 1 or a -1 shift


	# check 2nd frame for the +1 shift
	plus_1_stop = find_stop_codon(seq, drop_info[0]+1, len(seq))
	#print("frame two stop at position: " + str(plus_1_stop))
	# check 3rd frame for -1 shift
	minus_1_stop = find_stop_codon(seq, drop_info[0]-1, len(seq))
	#print("frame three stop at position: " + str(minus_1_stop))

	# get frame 1, 2 and frame 3 percentage before the drop
	percent_before_f1 = custom_get_percent_in_range(ss, 0, drop_info[0], "f1") 
	percent_before_f2 = custom_get_percent_in_range(ss, 0, drop_info[0], "f2") 
	percent_before_f3 = custom_get_percent_in_range(ss, 0, drop_info[0], "f3") 

	# get frame 2, 2 and frame 3 precentage after the drop
	percent_after_f1_f2 = custom_get_percent_in_range(ss,drop_info[0], plus_1_stop, "f1") 
	percent_after_f1_f3 = custom_get_percent_in_range(ss,drop_info[0], minus_1_stop, "f1") 
	percent_after_f2 = custom_get_percent_in_range(ss,drop_info[0], plus_1_stop, "f2") 
	percent_after_f3 = custom_get_percent_in_range(ss, drop_info[0], minus_1_stop, "f3") 
	#print("frame 2 percent before is: " + str(percent_before_f2) + " and after is " + str(percent_after_f2))
	#print("frame 3 percent before is: " + str(percent_before_f3) + " and after is " + str(percent_after_f3))
	#print("frame 1 precent before is: " + str(percent_before_f1) + " and after if frame 2 is " + str(percent_after_f1_f2) + " and after if frame 3 is " + str(percent_after_f1_f3))
	
	f2_increase = percent_after_f2 - percent_before_f2
	f3_increase = percent_after_f3 - percent_before_f3
	#print("f2 increase is " + str(f2_increase))
	#print("f3 increase is " + str(f3_increase))

	real_new_frame = "f2"
	real_increase = f2_increase
	real_length = plus_1_stop - drop_info[0]
	if (f3_increase > f2_increase):
		real_new_frame = "f3"
		real_increase = f3_increase	
		real_length = minus_1_stop - drop_info[0]

	if (real_length <  30):
		#outFile.write("Gene: " + current_gene + "\tNothing Found\n")
		continue

	before_percents = np.array(custom_get_range_of_percents(ss, 0, drop_info[0], "f1"))
	after_percents = np.array(custom_get_range_of_percents(ss, drop_info[0], drop_info[0]+real_length, "f1"))
	t_test = stats.ttest_ind(before_percents, after_percents, equal_var=False)
	p_value = t_test[1]
	if (p_value <= 0.0001):
		outFile.write("Gene: " + current_gene + "\tp-value: " + str(p_value) + "\tgene-length: " + str(len(seq)) + "\tshift-start-position: " + str(drop_info[0]) + "\tshift-length " + str(real_length) + "\n")
		f1_starts = "" 
		f_alt_starts = ""
		seq_info = ""
		for i in range(drop_info[0]-30,drop_info[0]+30,3):
			if (i == drop_info[0]):
				f1_starts += "| "
				f_alt_starts += "| "
				seq_info += "| "
			f1_starts += str(ss[i]) + "   "
			seq_info += seq[i:i+1] + " " + seq[i+1:i+2] + " " + seq[i+2:i+3] + " "
			if (real_new_frame == "f2"):
				f_alt_starts += " " + str(ss[i+1]) + "  " 
			else:
				f_alt_starts += "  " + str(ss[i-1]) + " "
		outFile.write(seq_info + "\n")
		outFile.write(f1_starts + "\n")
		outFile.write(f_alt_starts + "\n")
		percent_f1_all = ",".join(str(round(x,3)) for x in custom_get_range_of_percents(ss, 0, len(ss)-2, "f1"))
		outFile.write(percent_f1_all)
		outFile.write("\n")
		outFile.write("\n")
inFile.close()
outFile.close()
