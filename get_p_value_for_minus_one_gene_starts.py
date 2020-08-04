#!/usr/bin/python

import sys
import re
import random
# gene	stop	p_value	fraction_shifted	pattern_pvalue

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
for line in orfsFile:
	line = line.strip()
	if (line.startswith(">")):
		if (orf_name != ""):
			Orfs_Seq_Dict[orf_name] = dict() 
			Orfs_Seq_Dict[orf_name]["coding"] = orf_string
			if (orf_string[0:3] != "ATG"):
				print("Problem with header " + orf_string)
			orf_string = ""
		orf_name = line.split()[0][1:]
	else:
		orf_string += line
# process last orf
Orfs_Seq_Dict[orf_name] = dict() 
Orfs_Seq_Dict[orf_name]["coding"] = orf_string
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

genesWithMinusOneStartsDict = dict()
genesOutputInfo = dict()
minusDists = dict()
genesWithMinusOneStartsFile = open(sys.argv[5],'r')
geneMinusOneNames = []
all_dists = []
for line in genesWithMinusOneStartsFile:
	ts = line.strip().split(",")
	genesWithMinusOneStartsDict[ts[0]] = ts[4]
	genesOutputInfo[ts[0]] = ts
	geneMinusOneNames.append(ts[0])
genesWithMinusOneStartsFile.close()

genesOutFile = open(sys.argv[6], "w")


# now process the ids file
sgdToGeneName = dict()
inFile = open("saccharomyces_cerevisiae_gene_info.csv", "r")
inFile.readline()
for line in inFile:
    i_string = line.strip().split("\t")
    sgdToGeneName[i_string[1]] = i_string[3].split(" ")[0]
    

# now process the secondary structure file - this was the initial filter to see we can compare to PRFDB
inFile = open("saccharomyces_cerevisiae_mfe.csv", "r")
inFile.readline()
geneStructuresDict = dict()
for line in inFile:
	i_string = line.strip().split("\t")
	if (i_string[2] in sgdToGeneName):
		gene_name = sgdToGeneName[i_string[2]]
		start = int(i_string[4])
		if gene_name in geneStructuresDict:
			if start not in geneStructuresDict[gene_name]:
				geneStructuresDict[gene_name].append(start) 
		else:
			geneStructuresDict[gene_name] = []
			geneStructuresDict[gene_name].append(start) 
print(len(geneStructuresDict.keys()))
print(geneStructuresDict)


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

				for i in range(0,1):	
					# now compare to the AA sequence read in from the protein file
					if (Orfs_Seq_Dict[name]["strand"] == "-"):
						complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
						loc = random.randint(0,Orfs_Seq_Dict[name]["start"] - Orfs_Seq_Dict[name]["end"])
						if (name in genesWithMinusOneStartsDict and i == 0): # correct location for the real gene
							ts = genesOutputInfo[name]
							loc = int(ts[1])
					else: # positive
						loc = random.randint(0,Orfs_Seq_Dict[name]["end"] - Orfs_Seq_Dict[name]["start"])
						if (name in genesWithMinusOneStartsDict and i == 0): # correct location for the real gene
							ts = genesOutputInfo[name]
							loc = int(ts[1])
					secondary_structure_dist = 10000 
					if (name in geneStructuresDict):
						for idx in geneStructuresDict[name]:
							if (abs(loc - idx) < secondary_structure_dist):
								secondary_structure_dist = abs(loc - idx)

					all_dists.append(secondary_structure_dist)
					if (name in genesWithMinusOneStartsDict and i == 0): # correct location for the real gene
						minusDists[name] = secondary_structure_dist
						test_string = Ref_Dict[Orfs_Seq_Dict[name]["chrom"]][Orfs_Seq_Dict[name]["start"]:Orfs_Seq_Dict[name]["end"]]
						if (Orfs_Seq_Dict[name]["strand"] == "-"):
							test_string  = "".join(complement.get(base, base) for base in reversed(test_string)) # get the reverse complement
						slip = 100000 
						for test in ["AAA", "TTT"]:
							found_idxs = [m.start() for m in re.finditer(test, test_string)]
							for last_idx in found_idxs:
								if (abs(loc - last_idx) < slip):
									slip = last_idx	
						if (abs(loc - slip) < 30):
							if (slip%3 == 2):
								print("Name: " + name + " Slip: " + str(slip) + " found in correct frame\n")						
sorted_all_dists = sorted(all_dists)
print("len of all all_dists")
print(len(all_dists))


for name in geneMinusOneNames:
	dist = minusDists[name]	
	print(name)
	print(dist)
	idx = sorted_all_dists.index(dist)
	max_idx = max(loc for loc, val in enumerate(sorted_all_dists) if val == dist)
	p_value = random.randint(idx,max_idx)/len(all_dists)  
	print(str(dist) + " " + str(idx) + " " + str(max_idx) + " " +  str(p_value))
	ts = genesOutputInfo[name]
	if (p_value <= 0.05):  # THE P-value cutoff!!!!!
		genesOutFile.write(ts[0] + "," + ts[1] + "," + ts[2] + "," + ts[3] + "," + ts[4] + "," + str(round(float((ts[5])))) + "," + str(round(p_value,6)) + "\n")
print(len(sorted_all_dists))

my_set = set(sorted_all_dists)
new_list = list(my_set)
sorted_list = sorted(new_list)
print(sorted_list)
print(sorted_all_dists)
genesOutFile.close()

#TODO: add this file command and the ResultsForPaper.ipynb to READM and commit the files
pvaluesOutFile = open(sys.argv[7], "w")
for i in range(1, len(all_dists)):
	#p_value = i/len(all_sums)
	pvaluesOutFile.write(str(round(all_dists[i],6)) + ",")
pvaluesOutFile.write("\n")
pvaluesOutFile.close()
