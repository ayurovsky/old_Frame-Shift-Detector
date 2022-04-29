#!/usr/bin/python

import sys
import re
import random
# gene	stop	p_value	fraction_shifted	pattern_pvalue

random.seed(10)

###########################################################################
# First figure out what codons are
ctaa = dict() 
aatc = dict() 
codonsFile = open(sys.argv[1], 'r') 
for line in codonsFile:
	linesplit = line.split()
	ctaa[linesplit[0]] = linesplit[2]
	if (linesplit[2] in aatc):
		aatc[linesplit[2]].append(linesplit[0])
	else:
		aatc[linesplit[2]] = [linesplit[0]]

def aaseq(sequence):
	seq = "" 
	for i in range(0, len(sequence), 3):
		if len(sequence[i:i+3]) < 3:
			break
		seq += ctaa[sequence[i:i+3]]
	return seq

################################################################################
# create a map of ORF names
orfsFile = open(sys.argv[2],'r')
Orfs_Seq_Dict = dict()
orf_string = ""
orf_name = ""
orf_gene_name = ""
for line in orfsFile:
	line = line.strip()
	if (line.startswith(">")):
		if (orf_name != ""):
			Orfs_Seq_Dict[orf_name] = dict() 
			Orfs_Seq_Dict[orf_name]["coding"] = orf_string
			Orfs_Seq_Dict[orf_name]["gene_name"] = orf_gene_name
			if (orf_string[0:3] != "ATG"):
				print("Problem with header " + orf_string)
			orf_string = ""
		orf_name = line.split()[0][1:]
		orf_gene_name = line.split()[1]
	else:
		orf_string += line
# process last orf
Orfs_Seq_Dict[orf_name] = dict() 
Orfs_Seq_Dict[orf_name]["coding"] = orf_string
Orfs_Seq_Dict[orf_name]["gene_name"] = orf_gene_name
if (orf_string[0:3] != "ATG"):
	print("Problem with header " + orf_string)

print(str(len(Orfs_Seq_Dict)) + " protein coding orfs before any checks")


################################################################################
# read in the reference sequences for the chromosomes
refFile = open(sys.argv[3],'r')
ref_chr = ""
ref_seq = ""
Ref_Dict = dict()
for line in refFile:
	line = line.strip()
	if (line.startswith(">")):
		if (ref_chr != ""):
			Ref_Dict[ref_chr] = ref_seq.upper()
		ref_chr = line.split()[0][1:]
		ref_seq = ""	
	else:
		ref_seq += line
# last line	
Ref_Dict[ref_chr] = ref_seq.upper()

print("Read in the reference chromosomes\n")
################################################################################

# need to allow protein genes with programmed frameshifts:
#1. read in "gene" not "CDS";  read in the entire sequence, don't care if does not match the sequence in the file - read in the whole gene, use the stuff from the coordinates to match
#2. then if the same gene has introns, delete it for now...
#3. double-check that YIL009C-A(EST3) and YOR239W(ABP140) get in here correctly, and then make sure that they are not removed by intersections - remove other genes if necessary...

# read in the annotation file, check that protein AA sequence matches with the reference string codons, check that stop at the end
# introns no explicitly mentioned, filter out genes where ORF does not match the protein
annFile = open(sys.argv[4],'r')
Use_CDS_Dict = dict()
for line in annFile:
	if (line.startswith("#")):
		continue
	ts = line.strip().split()	
	name_search = re.compile(".*Name=([\(\)\.0-9a-zA-z_-]*);.*").search(ts[8])
	if (name_search != None):
		name = name_search.group(1)
		if (ts[2] == "five_prime_UTR_intron"):
			Use_CDS_Dict[name] = 1			

annFile.close()

genesWithPlusOneStartsDict = dict()
genesOutputInfo = dict()
plusDistance = dict()
plusSum = dict()
genesWithPlusOneStartsFile = open(sys.argv[5],'r')
genePlusOneNames = []
all_dists = []
all_sums = []
for line in genesWithPlusOneStartsFile:
	ts = line.strip().split(",")
	genesWithPlusOneStartsDict[ts[0]] = ts[4]
	genesOutputInfo[ts[0]] = ts
	genePlusOneNames.append(ts[0])
genesWithPlusOneStartsFile.close()

genesOutFile = open(sys.argv[6], "w")

found_gene = 0
annFile = open(sys.argv[4],'r')
for line in annFile:
	if (line.startswith("#")):
		continue
	ts = line.strip().split()	
	name_search = re.compile(".*Name=([\(\)\.0-9a-zA-z_-]*);.*").search(ts[8])
	if (name_search != None):
		name = name_search.group(1)
		if ((ts[2] == "gene" and name not in Use_CDS_Dict) or (ts[2] == "CDS" and name in Use_CDS_Dict)):
			if (name in Orfs_Seq_Dict): # only read in the genes that are in the full protein file, but do NOT check that the string is the same - this will allow the frame-shfited genes to pass through..
				# important - GFF is 1-based!!! the sam file read coordinates are also 1-based
				Orfs_Seq_Dict[name]["chrom"] = ts[0];
				Orfs_Seq_Dict[name]["start"] = int(ts[3])-1;
				Orfs_Seq_Dict[name]["end"] = int(ts[4])-1;
				Orfs_Seq_Dict[name]["strand"] = ts[6];
				if (Orfs_Seq_Dict[name]["strand"] == "-"):
					Orfs_Seq_Dict[name]["start"] = int(ts[4])-1; # for reverse direction gene, end is the start
					Orfs_Seq_Dict[name]["end"] = int(ts[3])-1; # for reverse direction gene, end is the start
				Orfs_Seq_Dict[name]["reads"] = 0;
				#Orfs_Seq_Dict[name]["parent_dbx_id"] = parent_gene_id = re.compile(".*,GeneID:([0-9]*),.*").search(ts[8]).group(1)

				for i in range(0,10):	
					# now compare to the AA sequence read in from the protein file
					test_string = ""
					limit = 60 
					if (Orfs_Seq_Dict[name]["strand"] == "-"):
						complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
						loc = random.randint(0,Orfs_Seq_Dict[name]["start"] - Orfs_Seq_Dict[name]["end"])
						if (name in genesWithPlusOneStartsDict and i == 0): # correct location for the real gene
							ts = genesOutputInfo[name]
							loc = int(ts[1])
						test_string = Ref_Dict[Orfs_Seq_Dict[name]["chrom"]][Orfs_Seq_Dict[name]["start"]-loc:Orfs_Seq_Dict[name]["start"]-loc+limit]
						test_string  = "".join(complement.get(base, base) for base in reversed(test_string)) # get the reverse complement
					else: # positive
						loc = random.randint(0,Orfs_Seq_Dict[name]["end"] - Orfs_Seq_Dict[name]["start"])
						if (name in genesWithPlusOneStartsDict and i == 0): # correct location for the real gene
							ts = genesOutputInfo[name]
							loc = int(ts[1])
							print("test left " + name)
						test_string = Ref_Dict[Orfs_Seq_Dict[name]["chrom"]][Orfs_Seq_Dict[name]["start"]+loc-limit:Orfs_Seq_Dict[name]["start"]+loc]
						if (name in genesWithPlusOneStartsDict and i == 0): # correct location for the real gene
							print("-" + test_string + "-")
							print(str(Orfs_Seq_Dict[name]["start"]+loc-limit) + " " + str(Orfs_Seq_Dict[name]["start"]+loc))
					oof_stop = limit 
					for test in ["TGA", "TAA", "TAG"]:
						found_idxs = [m.start() for m in re.finditer(test, test_string)]
						if (len(found_idxs)):
							last_idx = found_idxs[-1]
							if (limit - last_idx < oof_stop):
								oof_stop = limit - last_idx
					if (oof_stop > 3):
						all_dists.append(oof_stop)
					left_oof_stop = oof_stop
					# now try to the right of the random stop
					limit = 60
					if (Orfs_Seq_Dict[name]["strand"] == "-"):
						complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
						loc = random.randint(0,Orfs_Seq_Dict[name]["start"] - Orfs_Seq_Dict[name]["end"])
						if (name in genesWithPlusOneStartsDict and i == 0): # correct location for the real gene
							ts = genesOutputInfo[name]
							loc = int(ts[1])
						test_string = Ref_Dict[Orfs_Seq_Dict[name]["chrom"]][Orfs_Seq_Dict[name]["start"]-loc-limit:Orfs_Seq_Dict[name]["start"]-loc]
						test_string  = "".join(complement.get(base, base) for base in reversed(test_string)) # get the reverse complement
					else: # positive
						loc = random.randint(0,Orfs_Seq_Dict[name]["end"] - Orfs_Seq_Dict[name]["start"])
						if (name in genesWithPlusOneStartsDict and i == 0): # correct location for the real gene
							ts = genesOutputInfo[name]
							loc = int(ts[1])
							print("test right " + name)
						test_string = Ref_Dict[Orfs_Seq_Dict[name]["chrom"]][Orfs_Seq_Dict[name]["start"]+loc:Orfs_Seq_Dict[name]["start"]+loc+limit]
						if (name in genesWithPlusOneStartsDict and i == 0): # correct location for the real gene
							print("-" + test_string + "-")
							print(str(Orfs_Seq_Dict[name]["start"]+loc) + " " + str(Orfs_Seq_Dict[name]["start"]+loc+limit))
					oof_stop = limit 
					for test in ["TGA", "TAA", "TAG"]:
						found_idxs = [m.start() for m in re.finditer(test, test_string)]
						if (len(found_idxs)):
							first_idx = found_idxs[0]
							if (first_idx < oof_stop):
								oof_stop = first_idx
	
					#print("oof_stop " + str(oof_stop))
					if (name in genesWithPlusOneStartsDict and i == 0):
						plusDistance[name] = min(oof_stop, left_oof_stop)
						if (plusDistance[name] <= 3):
							all_dists.append(plusDistance[name])

				found_gene += 1
	#if (found_gene == 100):
	#	break

print(found_gene)
sorted_all_dists = sorted(all_dists)
print("len of all distances")
print(len(all_dists))
annFile.close()

codonSpeedsFile = open("codon_speeds.csv", 'r', encoding='utf-8-sig')
codonSpeedsDict = dict()
for line in codonSpeedsFile:
    ts = line.strip().split(",")
    codonSpeedsDict[ts[0]] = float(ts[1])
print(codonSpeedsDict)

found_gene = 0
annFile = open(sys.argv[4],'r')
for line in annFile:
	if (line.startswith("#")):
		continue
	ts = line.strip().split()	
	name_search = re.compile(".*Name=([\(\)\.0-9a-zA-z_-]*);.*").search(ts[8])
	if (name_search != None):
		name = name_search.group(1)
		if ((ts[2] == "gene" and name not in Use_CDS_Dict) or (ts[2] == "CDS" and name in Use_CDS_Dict)):
			if (name in Orfs_Seq_Dict): # only read in the genes that are in the full protein file, but do NOT check that the string is the same - this will allow the frame-shfited genes to pass through..
				# important - GFF is 1-based!!! the sam file read coordinates are also 1-based
				Orfs_Seq_Dict[name]["chrom"] = ts[0];
				Orfs_Seq_Dict[name]["start"] = int(ts[3])-1;
				Orfs_Seq_Dict[name]["end"] = int(ts[4])-1;
				Orfs_Seq_Dict[name]["strand"] = ts[6];
				if (Orfs_Seq_Dict[name]["strand"] == "-"):
					Orfs_Seq_Dict[name]["start"] = int(ts[4])-1; # for reverse direction gene, end is the start
					Orfs_Seq_Dict[name]["end"] = int(ts[3])-1; # for reverse direction gene, end is the start
				Orfs_Seq_Dict[name]["reads"] = 0;
				#Orfs_Seq_Dict[name]["parent_dbx_id"] = parent_gene_id = re.compile(".*,GeneID:([0-9]*),.*").search(ts[8]).group(1)

				for i in range(0,10):	
					# now compare to the AA sequence read in from the protein file
					test_string = ""
					limit = 60 
					if (Orfs_Seq_Dict[name]["strand"] == "-"):
						complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
						loc = random.randint(0,Orfs_Seq_Dict[name]["start"] - Orfs_Seq_Dict[name]["end"])
						if (name in genesWithPlusOneStartsDict and i == 0): # correct location for the real gene
							ts = genesOutputInfo[name]
							loc = int(ts[1])
						test_string = Ref_Dict[Orfs_Seq_Dict[name]["chrom"]][Orfs_Seq_Dict[name]["start"]-loc:Orfs_Seq_Dict[name]["start"]-loc+limit]
						test_string  = "".join(complement.get(base, base) for base in reversed(test_string)) # get the reverse complement
					else: # positive
						loc = random.randint(0,Orfs_Seq_Dict[name]["end"] - Orfs_Seq_Dict[name]["start"])
						if (name in genesWithPlusOneStartsDict and i == 0): # correct location for the real gene
							ts = genesOutputInfo[name]
							loc = int(ts[1])
							print("test left " + name)
						test_string = Ref_Dict[Orfs_Seq_Dict[name]["chrom"]][Orfs_Seq_Dict[name]["start"]+loc-limit:Orfs_Seq_Dict[name]["start"]+loc]
						if (name in genesWithPlusOneStartsDict and i == 0): # correct location for the real gene
							print("-" + test_string + "-")
							print(str(Orfs_Seq_Dict[name]["start"]+loc-limit) + " " + str(Orfs_Seq_Dict[name]["start"]+loc))
					speed_sum_left = 0
					for j in range(0, len(test_string), 3):
						codon = test_string[j:j+3]
						if codon in codonSpeedsDict:
							speed_sum_left += codonSpeedsDict[codon]

					all_sums.append(speed_sum_left)
					# now try to the right of the random stop
					limit = 60
					if (Orfs_Seq_Dict[name]["strand"] == "-"):
						complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
						loc = random.randint(0,Orfs_Seq_Dict[name]["start"] - Orfs_Seq_Dict[name]["end"])
						if (name in genesWithPlusOneStartsDict and i == 0): # correct location for the real gene
							ts = genesOutputInfo[name]
							loc = int(ts[1])
						test_string = Ref_Dict[Orfs_Seq_Dict[name]["chrom"]][Orfs_Seq_Dict[name]["start"]-loc-limit:Orfs_Seq_Dict[name]["start"]-loc]
						test_string  = "".join(complement.get(base, base) for base in reversed(test_string)) # get the reverse complement
					else: # positive
						loc = random.randint(0,Orfs_Seq_Dict[name]["end"] - Orfs_Seq_Dict[name]["start"])
						if (name in genesWithPlusOneStartsDict and i == 0): # correct location for the real gene
							ts = genesOutputInfo[name]
							loc = int(ts[1])
							print("test right " + name)
						test_string = Ref_Dict[Orfs_Seq_Dict[name]["chrom"]][Orfs_Seq_Dict[name]["start"]+loc:Orfs_Seq_Dict[name]["start"]+loc+limit]
						if (name in genesWithPlusOneStartsDict and i == 0): # correct location for the real gene
							print("-" + test_string + "-")
							print(str(Orfs_Seq_Dict[name]["start"]+loc) + " " + str(Orfs_Seq_Dict[name]["start"]+loc+limit))
	
					speed_sum_right = 0
					for j in range(0, len(test_string), 3):
						codon = test_string[j:j+3]
						if codon in codonSpeedsDict:
							speed_sum_right += codonSpeedsDict[codon]

					if (name in genesWithPlusOneStartsDict and i == 0):
						plusSum[name] = min(speed_sum_left,speed_sum_right)
						all_sums.append(plusSum[name])

				found_gene += 1
	#if (found_gene == 100):
	#	break

print(found_gene)
sorted_all_sums = sorted(all_sums)
print("len of all all_sums")
print(len(all_sums))


for name in genePlusOneNames:
	dist = plusDistance[name]	
	print(name)
	print(dist)
	idx = sorted_all_dists.index(dist)
	max_idx = max(loc for loc, val in enumerate(sorted_all_dists) if val == dist)
	p_value = random.randint(idx,max_idx)/len(all_dists) # dist 
	print(str(dist) + " " + str(idx) + " " + str(max_idx) + " " +  str(p_value))
	speed = plusSum[name]	
	idx = sorted_all_sums.index(speed)
	max_idx = max(loc for loc, val in enumerate(sorted_all_sums) if val == speed)
	speed_p_value = random.randint(idx,max_idx)/len(all_sums) # sum 
	print(str(speed) + " " + str(idx) + " " + str(max_idx) + " " +  str(speed_p_value))
	ts = genesOutputInfo[name]
#	if (p_value <= 0.05):  # THE P-value cutoff!!!!!
	genesOutFile.write(ts[0] + "," + Orfs_Seq_Dict[ts[0]]["gene_name"] + "," + ts[1] + "," + ts[2] + "," + ts[3] + "," + ts[4] + "," + str(round(float((ts[5])))) + "," + str(round(p_value,6)) + "," + str(round(speed_p_value,6)) + "\n")
print(len(sorted_all_dists))

my_set = set(sorted_all_dists)
new_list = list(my_set)
sorted_list = sorted(new_list)
print(sorted_list)
print(sorted_all_dists)
genesOutFile.close()
