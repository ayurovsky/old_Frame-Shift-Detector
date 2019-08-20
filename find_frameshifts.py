#!/usr/bin/python
import sys
import random
import scipy
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

inFile = open(sys.argv[1],'r')
gene = inFile.readline().strip()
start_idx = int(sys.argv[2])
end_idx = int(sys.argv[3])
outFile = open("output/" + sys.argv[2] + "_" + sys.argv[3] + "_" + sys.argv[1], 'w')

counter = 0
while gene:
	seq = inFile.readline()
	starts_string = inFile.readline()
	current_gene = gene
	gene = inFile.readline().strip()
#	if (not (current_gene == "YIL009C-A" or current_gene == "YOR239W")):
#		continue
	if (counter < start_idx or counter >= end_idx):
		counter += 1
		continue
	else:
		counter += 1

	ss = starts_string.strip()[:-1].split(",")

	# get the max drop - % change in frame 1 reads	
	#print("_" + current_gene + "_")
	#print(seq)
	drop_info = get_drop(seq, ss)
	#print("Max drop location for running 11 avearge of percent location is: " + str(drop_info[0]) + " " + str(drop_info[1]))

	# decide if this is a + 1 or a -1 shift


	# check 2nd frame for the +1 shift
	plus_1_stop = find_stop_codon(seq, drop_info[0]+1, len(seq))
	#print("frame two stop at position: " + str(plus_1_stop))
	# check 3rd frame for -1 shift
	minus_1_stop = find_stop_codon(seq, drop_info[0]-1, len(seq))
	#print("frame three stop at position: " + str(minus_1_stop))

	# get frame 2 and frame 3 percentage before the drop
	percent_before_f1 = custom_get_percent_in_range(ss, 0, drop_info[0], "f1") 
	percent_before_f2 = custom_get_percent_in_range(ss, 0, drop_info[0], "f2") 
	percent_before_f3 = custom_get_percent_in_range(ss, 0, drop_info[0], "f3") 

	# get frame 2 and frame 3 precentage after the drop
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

	real_new_frame = "f2"	# TODO: check for -1s
	real_increase = f2_increase
	real_length = plus_1_stop - drop_info[0]
	if (f3_increase > f2_increase):
		real_new_frame = "f3"
		real_increase = f3_increase	
		real_length = minus_1_stop - drop_info[0]

	if (real_length <  30):
		#outFile.write("Gene: " + current_gene + "\tNothing Found\n")
		continue
	#print("New frame is " + real_new_frame + " and new frame increase is " + str(real_increase) + " and length of new tail is " + str(real_length))
	#print("Average before drop : " + str(mean(drop_info[2][:int(drop_info[0]/3)])))
	#print("Average after drop : " + str(mean(drop_info[2][int(drop_info[0]/3):int(drop_info[0]/3+real_length/3)])))
	real_difference_in_averages = mean(drop_info[2][:int(drop_info[0]/3)]) - mean(drop_info[2][int(drop_info[0]/3):int(drop_info[0]/3+real_length/3)])

	num_permutations = 1000
	num_better = 0 	
	for i in range(0,num_permutations):
		scramble_starts = ss[:]
		random.shuffle(scramble_starts)		
		scramble_drop_info = get_drop(seq, scramble_starts)
		
		plus_1_stop = find_stop_codon(seq, scramble_drop_info[0]+1, len(seq)) #TODO: checck for -1s
		minus_1_stop = find_stop_codon(seq, scramble_drop_info[0]-1, len(seq))
		percent_before_f2 = custom_get_percent_in_range(ss, 0, scramble_drop_info[0], "f2") 
		percent_before_f3 = custom_get_percent_in_range(ss, 0, scramble_drop_info[0], "f3") 
		percent_after_f2 = custom_get_percent_in_range(ss, scramble_drop_info[0], plus_1_stop, "f2") 
		percent_after_f3 = custom_get_percent_in_range(ss, scramble_drop_info[0], minus_1_stop, "f3") 
		f2_increase = percent_after_f2 - percent_before_f2
		f3_increase = percent_after_f3 - percent_before_f3
		new_frame = "f2"	
		increase = f2_increase
		length = plus_1_stop - scramble_drop_info[0]
		if (f3_increase > f2_increase):
			new_frame = "f3"
			increase = f3_increase	
			length = minus_1_stop - scramble_drop_info[0]
		if (length >= 30):
			#print("Scramble frame is " + new_frame + " and new frame increase is " + str(increase) + " and length of new tail is " + str(length))
			#print("Average before drop : " + str(mean(scramble_drop_info[2][:int(scramble_drop_info[0]/3)])) + "\n")
			#print("Average after drop : " + str(mean(scramble_drop_info[2][int(scramble_drop_info[0]/3):int(scramble_drop_info[0]/3+length/3)])) + "\n")
			difference_in_averages = mean(scramble_drop_info[2][:int(scramble_drop_info[0]/3)]) - mean(scramble_drop_info[2][int(scramble_drop_info[0]/3):int(scramble_drop_info[0]/3+length/3)])
			if (difference_in_averages >= real_difference_in_averages):
				num_better += 1 
	p_value = -1 
	if (num_better == 0):
		p_value = 1.0/num_permutations
	else:
		p_value = num_better/num_permutations
	#outFile.write("Gene: " + current_gene + "\tp-value: " + str(p_value) + "\n")
	if (p_value <= 0.005):
		outFile.write("Gene: " + current_gene + "\tp-value: " + str(p_value) + "\tgene-length: " + str(len(seq)) + "\tshift-start-position: " + str(drop_info[0]) + "\tshift-length " + str(real_length) + "\n")
		f1_starts = "" 
		f_alt_starts = ""
		seq_info = ""
		for i in range(drop_info[0]-30,drop_info[0]+30,3):
			f1_starts += str(ss[i]) + "\t"
			seq_info += seq[i:i+3] + "\t"
			if (real_new_frame == "f2"):
				f_alt_starts += " " + str(ss[i+1]) + "\t"
			else:
				f_alt_starts += "  " + str(ss[i-1]) + "\t"
			if (i == drop_info[0]):
				f1_starts += "|\t"
				f_alt_starts += "|\t"
				seq_info += "|\t"
		outFile.write(seq_info + "\n")
		outFile.write(f1_starts + "\n")
		outFile.write(f_alt_starts + "\n")
		outFile.write("\n")
inFile.close()
outFile.close()
