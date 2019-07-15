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

################################################################################

# read in the annotation file, check that protein AA sequence matches with the reference string codons, check that stop at the end
# introns no explicitly mentioned, filter out genes where ORF does not match the protein
annFile = open(sys.argv[4],'r')
Exclude_Dict = dict()
for line in annFile:
	if (line.startswith("#")):
		continue
	ts = line.strip().split()	
	if (ts[2] == "CDS" and "Name=" in ts[8]):
		name  = re.compile(".*Name=([\(\)\.0-9a-zA-z_-]*);.*").search(ts[8]).group(1)
		if (name in Orfs_Seq_Dict):
			# important - GFF is 1-based!!! the sam file read coordinates are also 1-based
			Orfs_Seq_Dict[name]["chrom"] = ts[0];
			Orfs_Seq_Dict[name]["start"] = int(ts[3])-1;
			Orfs_Seq_Dict[name]["end"] = int(ts[4])-1;
			Orfs_Seq_Dict[name]["strand"] = ts[6];
			if (Orfs_Seq_Dict[name]["strand"] == "-"):
				Orfs_Seq_Dict[name]["start"] = int(ts[4])-1; # for reverse direction gene, end is the start
				Orfs_Seq_Dict[name]["end"] = int(ts[3])-1; # for reverse direction gene, end is the start
			Orfs_Seq_Dict[name]["reads"] = 0;
			Orfs_Seq_Dict[name]["parent_dbx_id"] = parent_gene_id = re.compile(".*,GeneID:([0-9]*),.*").search(ts[8]).group(1)

			# now compare to the AA sequence read in from the protein file
			full_string = ""
			if (Orfs_Seq_Dict[name]["strand"] == "-"):
				full_string = Ref_Dict[Orfs_Seq_Dict[name]["chrom"]][Orfs_Seq_Dict[name]["end"]:Orfs_Seq_Dict[name]["start"]+1]
				complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
				full_string  = "".join(complement.get(base, base) for base in reversed(full_string)) # get the reverse complement
			else:
				full_string = Ref_Dict[Orfs_Seq_Dict[name]["chrom"]][Orfs_Seq_Dict[name]["start"]:Orfs_Seq_Dict[name]["end"]+1]
			if (Orfs_Seq_Dict[name]["aa"] != aaseq(full_string)[:-1]): # check that annotated proteins match the reference sequence translated to aa
				#print("problem with " + name + " annotated protein sequence not matching with reference sequence translated to aa ->  introns?")
				#print(aaseq(full_string))
				#print(Orfs_Seq_Dict[name]["aa"])
				Orfs_Seq_Dict.pop(name)
			elif (aaseq(full_string)[-1:] != "*"):
				print("problem with " + name + " protein not ending with a stop codon!")
				Orfs_Seq_Dict.pop(name)
			else: # save the seq?
				Orfs_Seq_Dict[name]["seq"] = full_string
	# this is a top-level structure - will catch all genes for proteins, tRNAs, rRNAs, but not regions
	if ("Parent=" not in ts[8]):
		i_d  = re.compile(".*ID=([\(\)\.0-9a-zA-z:_-]*);.*").search(ts[8]).group(1)
		if (":" not in i_d):
			dbx_id = re.compile(".*GeneID:([0-9]*).*").search(ts[8]).group(1)
			Exclude_Dict[dbx_id] = ts

print(str(len(Orfs_Seq_Dict)) + " protein coding orfs after checking for introns and other mis-matches")
################################################################################

# now remove proteins that intersect with other proteins or tRNA, rRNA
for orf_name in list(Orfs_Seq_Dict):
	if "parent_dbx_id" in Orfs_Seq_Dict[orf_name]: # orf has been found, parent dbx_id identified
		for toplevel_dbx_id in Exclude_Dict:
			if toplevel_dbx_id != Orfs_Seq_Dict[orf_name]["parent_dbx_id"]: # not self TODO: maybe here get the real gene name, if useful for later..
				orf_start = min(Orfs_Seq_Dict[orf_name]["start"], Orfs_Seq_Dict[orf_name]["end"])
				orf_end = max(Orfs_Seq_Dict[orf_name]["start"], Orfs_Seq_Dict[orf_name]["end"])
				top_start = int(Exclude_Dict[toplevel_dbx_id][3])
				top_end = int(Exclude_Dict[toplevel_dbx_id][4])
				top_chr = Exclude_Dict[toplevel_dbx_id][0]
				if (Orfs_Seq_Dict[orf_name]["chrom"] == top_chr): 
					if (((orf_end >= top_start) and (orf_end <= top_end)) or ((orf_start >= top_start) and (orf_start <= top_end)) or ((orf_start <= top_start) and (orf_end >= top_end)) or ((orf_start >= top_start) and (orf_end <= top_end))):  #overlap
						#print("deleting ") 
						#print(Orfs_Seq_Dict[orf_name])
						#print(str(orf_start) + " " + str(orf_end))
						#print("intersected with ")
						#print(Exclude_Dict[toplevel_dbx_id])
						#print(str(top_start) + " " + str(top_end))
						del Orfs_Seq_Dict[orf_name]
						break
	else:
		del Orfs_Seq_Dict[orf_name]

print(str(len(Orfs_Seq_Dict)) + " protein coding orfs after intersect removals")

################################################################################
# assign every position in an orf to its orf name
Pos_Orf_Dict = dict()
for orf_name in Orfs_Seq_Dict:
	orf_start = min(Orfs_Seq_Dict[orf_name]["start"], Orfs_Seq_Dict[orf_name]["end"])
	orf_end = max(Orfs_Seq_Dict[orf_name]["start"], Orfs_Seq_Dict[orf_name]["end"])
	chrom = Orfs_Seq_Dict[orf_name]["chrom"]
	if (chrom not in Pos_Orf_Dict):
		Pos_Orf_Dict[chrom] = dict() 
	for pos in range(orf_start,orf_end+1):
		Pos_Orf_Dict[chrom][pos] = orf_name


################################################################################
# now process the sam file
