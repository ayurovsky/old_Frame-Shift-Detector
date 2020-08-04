#!/usr/bin/python

import sys
import re

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

genesWithBeforeGeneStartsDict = dict()
genesOutputInfo = dict()
beforeStartsDistance = dict()
genesWithBeforeGeneStartsFile = open(sys.argv[5],'r')
geneBeforeStartNames = []
for line in genesWithBeforeGeneStartsFile:
	ts = line.strip().split(",")
	genesWithBeforeGeneStartsDict[ts[0]] = ts[4]
	genesOutputInfo[ts[0]] = ts
	geneBeforeStartNames.append(ts[0])
genesWithBeforeGeneStartsFile.close()

genesOutFile = open(sys.argv[6], "w")

all_starts = []
before_gene_starts = []
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

				# now compare to the AA sequence read in from the protein file
				test_string = ""
				#full_string = ""
				limit = 360 
				if (Orfs_Seq_Dict[name]["strand"] == "-"):
					complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
					#full_string = Ref_Dict[Orfs_Seq_Dict[name]["chrom"]][Orfs_Seq_Dict[name]["end"]:Orfs_Seq_Dict[name]["start"]+1]
					#full_string  = "".join(complement.get(base, base) for base in reversed(full_string)) # get the reverse complement
					test_string = Ref_Dict[Orfs_Seq_Dict[name]["chrom"]][Orfs_Seq_Dict[name]["start"]:Orfs_Seq_Dict[name]["start"]+limit]
					test_string  = "".join(complement.get(base, base) for base in reversed(test_string)) # get the reverse complement
				else: # positive
					#full_string = Ref_Dict[Orfs_Seq_Dict[name]["chrom"]][Orfs_Seq_Dict[name]["start"]:Orfs_Seq_Dict[name]["end"]+1]
					test_string = Ref_Dict[Orfs_Seq_Dict[name]["chrom"]][Orfs_Seq_Dict[name]["start"]-limit:Orfs_Seq_Dict[name]["start"]]
				for_plus_one = limit 
				for_minus_one = limit 
				print(test_string)
				print([m.start() for m in re.finditer("ATG", test_string)])					
				print([m.start() for m in re.finditer("TTG", test_string)])					
				print([m.start() for m in re.finditer("ATA", test_string)])					
				for test in ["ATG", "TTG", "ATA"]:
					found_idxs = [m.start() for m in re.finditer(test, test_string)]
					if (len(found_idxs)):
						last_idx = found_idxs[-1]
						if (last_idx % 3 == 1): # for plus_one
							if (limit - last_idx < for_plus_one):
								for_plus_one = limit - last_idx
						elif (last_idx % 3 == 2): # for minus_one
							if (limit - last_idx < for_minus_one):
								for_minus_one = limit - last_idx
				print("plus one " + str(for_plus_one))
				print("minus one " + str(for_minus_one))
				all_starts.append(for_plus_one)
				all_starts.append(for_minus_one)	
				if (name in genesWithBeforeGeneStartsDict):
					if genesWithBeforeGeneStartsDict[name] == "plus":
						before_gene_starts.append(for_plus_one)
						beforeStartsDistance[name] = for_plus_one
					else:
						before_gene_starts.append(for_minus_one)
						beforeStartsDistance[name] = for_minus_one

				found_gene += 1
	#if (found_gene == 100):
	#	break

print(found_gene)
sorted_all_starts = sorted(all_starts)

for name in geneBeforeStartNames:
	start = beforeStartsDistance[name]	
	idx = max(1, (sorted_all_starts.index(start)))
	p_value = start/len(all_starts)
	print(str(start) + " " + str(idx) + " " + str(p_value))
	ts = genesOutputInfo[name]
	genesOutFile.write(ts[0] + "," + ts[1] + "," + ts[2] + "," + ts[3] + "," + str(round(p_value,6)) + "\n")
print(len(sorted_all_starts))
print(len(before_gene_starts))
print(before_gene_starts)

genesOutFile.close()
pvaluesOutFile = open(sys.argv[7], "w")
for i in range(1, len(all_starts)):
    p_value = i/len(all_starts)
    pvaluesOutFile.write(str(round(p_value,6)) + ",")
pvaluesOutFile.write("\n") 
pvaluesOutFile.close()
#TODO: add this file command and the ResultsForPaper.ipynb to READM and commit the files

