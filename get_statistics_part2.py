#!/usr/bin/python
import sys

def maximum(a, b, c): 
	if (a >= b) and (a >= c): 
		largest = a 
	elif (b >= a) and (b >= c): 
		largest = b 
	else: 
		largest = c 
	return largest 

inFile = open(sys.argv[1],'r')
gene = inFile.readline()
num_windows = 0
# all bigger_stats code from first version
#f1_bigger_stats = "" 
#f2_bigger_stats = "" 
#f3_bigger_stats = "" 

##max_gene_name = ""
##max_gene_rpkm = 0 
f1_higher_dict = dict()
f2_higher_dict = dict()
f3_higher_dict = dict()
total = 0
while gene:
	seq = inFile.readline()
	starts_string = inFile.readline()
	ss = starts_string.split(",")
	gene_total_reads = 0
	
	for i in range(0, len(seq) - 3, 3): # i is in frame
		f1 = int(ss[i])
		f2 = int(ss[i+1])
		f3 = int(ss[i+2])
		all_f = f1 + f2 + f3
		if (all_f ==0):
			continue
		total += all_f
		bucket = int(all_f/5)
		if (all_f >= 100):
			bucket = 20
		if (f1 == maximum(f1, f2, f3)):	
			if bucket not in f1_higher_dict:
				f1_higher_dict[bucket] = 0 
			f1_higher_dict[bucket] += 1
		elif (f2 == maximum(f1, f2, f3)):	
			if bucket not in f2_higher_dict:
				f2_higher_dict[bucket] = 0 
			f2_higher_dict[bucket] += 1
		else:	
			if bucket not in f3_higher_dict:
				f3_higher_dict[bucket] = 0 
			f3_higher_dict[bucket] += 1
	gene = inFile.readline()
	##for i in range(0, len(seq) - 3, 3): # i is in frame
	##gene_total_reads += int(ss[i]) + int(ss[i+1]) + int(ss[i+2])
	##	gene_rpkm = float(gene_total_reads)/len(seq)
	##if (gene_rpkm > max_gene_rpkm):
	##	max_gene_rpkm = gene_rpkm
	##	max_gene_name = gene 
	#f1_bigger = 0
	#f2_bigger = 0
	#f3_bigger = 0
	#for i in range(0, len(seq) - 30, 3): # i is in frame
		# all bigger_stats code from first version
		# check that there are at least 60 reads in a 10 codon window
		#starts_in_window = 0
		#for j in range(0,30):
		#	starts_in_window += int(ss[i+j])
		#if (starts_in_window >= 120):
			# printing only to get the statistics for a graph of sliding windows
			#print_string = ""
			#for j in range(0,30):      
			#	print_string += str(ss[i+j]) + ","
			#print(print_string)
		#	num_windows += 1
			# now check how the periodicity is going - get a single codon in the middle, and the check the reads
		#	f1_starts = int(ss[i+15])
		#	f2_starts = int(ss[i+16])
		#	f3_starts = int(ss[i+17])
		#	if (f1_starts == maximum(f1_starts, f2_starts, f3_starts)):
		#		f1_bigger += 1
		#	elif (f2_starts == maximum(f1_starts, f2_starts, f3_starts)):
		#		f2_bigger += 1
		#	else:
		#		f3_bigger += 1
	
	#print(str(f1_bigger) + " " + str(f2_bigger) + " " + str(f3_bigger))
	# all bigger_stats code from first version
	#if ((f1_bigger + f2_bigger + f3_bigger) > 10):
	#	f1_bigger_stats += str(round(float(f1_bigger)/(f1_bigger + f2_bigger + f3_bigger),3)) + ","
	#	f2_bigger_stats += str(round(float(f2_bigger)/(f1_bigger + f2_bigger + f3_bigger),3)) + ","
	#	f3_bigger_stats += str(round(float(f3_bigger)/(f1_bigger + f2_bigger + f3_bigger),3)) + ","
#print(f1_bigger_stats)
#print(f2_bigger_stats)
#print(f3_bigger_stats)
##print(max_gene_name)
##print(round(max_gene_rpkm))
print(total)
f1_string = ""
for i in range(0,21):
	f1_string += str(round(f1_higher_dict[i]/(f1_higher_dict[i] + f2_higher_dict[i] + f3_higher_dict[i]), 3)) + ","
print(f1_string)
f2_string = ""
for i in range(0,21):
	f2_string += str(round(f2_higher_dict[i]/(f1_higher_dict[i] + f2_higher_dict[i] + f3_higher_dict[i]), 3)) + ","
print(f2_string)
f3_string = ""
for i in range(0,21):
	f3_string += str(round(f3_higher_dict[i]/(f1_higher_dict[i] + f2_higher_dict[i] + f3_higher_dict[i]), 3)) + ","
print(f3_string)
inFile.close()
