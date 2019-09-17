#!/usr/bin/python
import sys
import random

def find_stop_codon(sequence, start, end):
	stop = -1
	for i in range(start,end,3):
		if (sequence[i:i+3] in ["TGA", "TAA", "TAG"]):
			stop = i
			break
	return stop

# end of functions

inFile = open(sys.argv[1],'r')
gene = inFile.readline().strip()

outFile = open(sys.argv[1] + "_with_simulated_shifts", 'w')
infoFile = open(sys.argv[1] + "_simulated_shifts_info", 'w')

simulated_shifts_dict = dict()
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
	# select the prob for the shift - this is the probability with which the ribosome shifts
	shift_prob = int(random.random()*10)/10
	if (shift_prob == 0.0):
		shift_prob = 0.05
#	print(shift_prob)
	if (shift_prob not in simulated_shifts_dict):
		simulated_shifts_dict[shift_prob] = 0 
	simulated_shifts_dict[shift_prob] += 1 

	if (len(ss)%3 != 0):
		#TODO!!!! in the get_statistics step, fix the genes with the 5_PRIME_UTR - remove the sequence in the beginning!!!!
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

	print(str(shift_direction) + " " + str(shift_location) + " " + str(stop_location))

	# assume the ribosome "slips", so all the read that were previously in frame will be shifted
 	# assume the width of the ribosome footprint is always 30, so the "visible start" of the shift is shift_location-15 and the visible end is stop_location-15: TODO: make more messy/probabilistic later		
	new_ss = [0] * len(ss)
	functional_start = shift_location - 15
	functional_end = stop_location - shift_direction - 15 # to get into frame 1 again
	# first fill in the reads before the shift
	for i in range(0, functional_start):
		new_ss[i] = ss[i]
	# now fill in the reads in the shift interval - they move with shift_prob, and stay otherwise
	for i in range(functional_start, functional_end):
		for j in range(0,ss[i]):
			if (random.random() < shift_prob):
				new_ss[i+shift_direction] += 1
			else:
				new_ss[i] += 1
	# now fill in reads after the shift - should be 1 - shift_prob for each read
	for i in range(functional_end,len(ss)):
		# for each read, keep it with given probability
		for j in range(0,ss[i]):
			if (random.random() < (1.0 - shift_prob)):
				new_ss[i] += 1
		
	
	#TODO: assuming all above is correct, write out the new_ss, the shift information, and display
	outFile.write(",".join([str(x) for x in new_ss]) + ",\n")
	infoFile.write(current_gene + "," + str(shift_location) + "," + str(stop_location) + "," + str(shift_direction) + "," + str(shift_prob) + "\n")
outFile.close()
