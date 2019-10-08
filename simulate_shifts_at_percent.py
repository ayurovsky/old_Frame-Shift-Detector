#!/usr/bin/python
import sys
import random
import math

def find_stop_codon(sequence, start, end):
	stop = -1
	for i in range(start,end,3):
		if (sequence[i:i+3] in ["TGA", "TAA", "TAG"]):
			stop = i
			break
	return stop

# end of functions

inFile = open(sys.argv[1],'r')
shift_prob = float(sys.argv[2])
gene = inFile.readline().strip()

outFile = open(sys.argv[1] + "_with_simulated_shifts_at_" + sys.argv[2], 'w')
infoFile = open(sys.argv[1] + "_simulated_shifts_info_at_" + sys.argv[2], 'w')

while gene:
	seq = inFile.readline()
	starts_string = inFile.readline()
	current_gene = gene
	gene = inFile.readline().strip()
	ss = [int(x) for x in starts_string.strip()[:-1].split(",")]
	rpbm = sum(ss)/len(ss)
	# output the first two lines - don't change
	outFile.write(current_gene + "\n")
	outFile.write(seq)
	if (rpbm < 10.0 or (current_gene == "YIL009C-A" or current_gene == "YOR239W")):
		outFile.write(",".join([str(x) for x in ss]) + ",\n")
		continue
	# ~600 genes > 10.0
	#print(shift_prob)

	if (len(ss)%3 != 0):
		#print("problem with gene " + current_gene)
		outFile.write(",".join([str(x) for x in ss]) + ",\n")
		continue 

	shift_location = -1
	shift_direction = 0
	stop_location = -1	
	try_counter = 0
	while ((stop_location - shift_location < 50) and (try_counter < 100)):
		try_counter += 1
		# now randomly select the location of the shift
		shift_location = random.randrange(11,(len(ss)-30)/3)*3
		# randomly select shift direction
		shift_direction = 1
		if (random.random() < 0.5):
			shift_direction = -1
		# get the stop codon in that frame
		stop_location = find_stop_codon(seq, shift_location + shift_direction, len(seq)) 
	if (try_counter == 100): # give up on this gene
		outFile.write(",".join([str(x) for x in ss]) + ",\n")
		continue

	#print(str(shift_direction) + " " + str(shift_location) + " " + str(stop_location))

	# assume the ribosome "slips", so all the read that were previously in frame will be shifted
 	# assume the width of the ribosome footprint is always 30, so the "visible start" of the shift is shift_location-5 and the visible end is stop_location-5: TODO: make more messy/probabilistic later		
	new_ss = [0] * len(ss)
	functional_start = shift_location - 5
	functional_end = stop_location - shift_direction - 5 # to get into frame 1 again
	# first fill in all the reads 
	for i in range(0, len(ss)):
		new_ss[i] = ss[i]
	# now fill in the reads in the shift interval - they get augmented with the shift probability 
	for i in range(functional_start, functional_end):
		#print(str(ss[i]) + " " + str(new_ss[i+shift_direction]))
		for j in range(0,ss[i]):
			if (random.random() < (shift_prob%1.0)): # augment with probability < 1.0
				 new_ss[i+shift_direction] += 1
			new_ss[i+shift_direction] += math.floor(shift_prob) # augment for probabilty >= 1 : ex. 100 percent, 200 percent
		#print(str(ss[i]) + " " + str(new_ss[i+shift_direction]))
	#print(ss[functional_start:functional_end])	
	#print(new_ss[functional_start:functional_end])	
	#TODO: assuming all above is correct, write out the new_ss, the shift information, and display
	outFile.write(",".join([str(x) for x in new_ss]) + ",\n")
	infoFile.write(current_gene + "," + str(shift_location) + "," + str(stop_location) + "," + str(shift_direction) + "," + str(shift_prob) + "\n")
outFile.close()
