#!/usr/bin/python

import sys

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
aa_string = ""
orf_name = ""
for line in orfsFile:
	line = line.strip()
	if (line.startswith(">")):
		if (orf_name != ""):
			Orfs_Seq_Dict[orf_name] = dict() 
			Orfs_Seq_Dict[orf_name]["aa"] = aa_string
			if (aa_string[0] != "M"):
				print("Problem with header " + aa_string)
			aa_string = ""
		orf_name = line.split()[0][1:]
	else:
		aa_string += line
# process last orf
Orfs_Seq_Dict[orf_name] = dict() 
Orfs_Seq_Dict[orf_name]["aa"] = aa_string
if (aa_string[0] != "M"):
	print("Problem with header " + aa_string)


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

################################################################################

# read in the annotation file, remove genes with introns, check that protein AA sequence matches with the reference string codons, check that stop at the end
annFile = open(sys.argv[4],'r')
for line in annFile:
	fields = line.strip().split()	


