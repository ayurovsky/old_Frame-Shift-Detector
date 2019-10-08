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
annFile = open(sys.argv[4],'r')

Exclude_Dict = dict()
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
				full_string = ""
				if (Orfs_Seq_Dict[name]["strand"] == "-"):
					full_string = Ref_Dict[Orfs_Seq_Dict[name]["chrom"]][Orfs_Seq_Dict[name]["end"]:Orfs_Seq_Dict[name]["start"]+1]
					complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
					full_string  = "".join(complement.get(base, base) for base in reversed(full_string)) # get the reverse complement
				else:
					full_string = Ref_Dict[Orfs_Seq_Dict[name]["chrom"]][Orfs_Seq_Dict[name]["start"]:Orfs_Seq_Dict[name]["end"]+1]
				if (name == "YIL009C-A" or name == "YOR239W"):
					print("for " + name + " protein seq is " + str(len(Orfs_Seq_Dict[name]["coding"])) + " while genomic seq is " + str(len(full_string)) + "\n")
				Orfs_Seq_Dict[name]["seq"] = full_string
		elif (ts[2] == "intron"):
			if (name in Orfs_Seq_Dict):
				Orfs_Seq_Dict.pop(name)
		# this is a top-level structure - will catch all genes for proteins, tRNAs, rRNAs, but not regions
		if ("Parent=" not in ts[8]):
			Exclude_Dict[name] = ts
			#i_d  = re.compile(".*ID=([\(\)\.0-9a-zA-z:_-]*);.*").search(ts[8]).group(1)
			#if (":" not in i_d):
			#	dbx_id = re.compile(".*GeneID:([0-9]*).*").search(ts[8]).group(1)
		#	Exclude_Dict[dbx_id] = ts

print(str(len(Orfs_Seq_Dict)) + " protein coding orfs after checking for introns and other mis-matches")

if ("YIL009C-A" not in Orfs_Seq_Dict):
	print("YIL009C-A not in Orfs_Seq_Dict")
if ("YOR239W" not in Orfs_Seq_Dict):
	print("YOR239W not in Orfs_Seq_Dict")
################################################################################

# now remove proteins that intersect with other proteins or tRNA, rRNA
for orf_name in list(Orfs_Seq_Dict):
	if "chrom" in Orfs_Seq_Dict[orf_name]: # orf has been found, parent dbx_id identified
		for toplevel_struct in Exclude_Dict:
			if toplevel_struct != orf_name: # not self TODO: maybe here get the real gene name, if useful for later..
				orf_start = min(Orfs_Seq_Dict[orf_name]["start"], Orfs_Seq_Dict[orf_name]["end"])
				orf_end = max(Orfs_Seq_Dict[orf_name]["start"], Orfs_Seq_Dict[orf_name]["end"])
				top_start = int(Exclude_Dict[toplevel_struct][3])
				top_end = int(Exclude_Dict[toplevel_struct][4])
				top_chr = Exclude_Dict[toplevel_struct][0]
				if (Orfs_Seq_Dict[orf_name]["chrom"] == top_chr): 
					if (((orf_end >= top_start) and (orf_end <= top_end)) or ((orf_start >= top_start) and (orf_start <= top_end)) or ((orf_start <= top_start) and (orf_end >= top_end)) or ((orf_start >= top_start) and (orf_end <= top_end))):  #overlap
						if (orf_name == "YIL009C-A" or orf_name == "YOR239W"):
							print("deleting " + orf_name + " because it intersects with " + toplevel_struct)
						#print("deleting ") 
						#print(Orfs_Seq_Dict[orf_name])
						#print(str(orf_start) + " " + str(orf_end))
						#print("intersected with ")
						#print(Exclude_Dict[toplevel_struct])
						#print(str(top_start) + " " + str(top_end))
						del Orfs_Seq_Dict[orf_name]
						break
	else:
		del Orfs_Seq_Dict[orf_name]

print(str(len(Orfs_Seq_Dict)) + " protein coding orfs after intersect removals")

if ("YIL009C-A" not in Orfs_Seq_Dict):
	print("YIL009C-A not in Orfs_Seq_Dict after intersect removals")
if ("YOR239W" not in Orfs_Seq_Dict):
	print("YOR239W not in Orfs_Seq_Dict after intersect removals")

################################################################################
# assign every position in an orf to its orf name
Pos_Orf_Dict = dict()
for orf_name in Orfs_Seq_Dict:
	if (orf_name == "YIL009C-A" or orf_name == "YOR239W"):
		print("assigning pos orf to " + orf_name)
	orf_start = min(Orfs_Seq_Dict[orf_name]["start"], Orfs_Seq_Dict[orf_name]["end"])
	orf_end = max(Orfs_Seq_Dict[orf_name]["start"], Orfs_Seq_Dict[orf_name]["end"])
	chrom = Orfs_Seq_Dict[orf_name]["chrom"]
	if (chrom not in Pos_Orf_Dict):
		Pos_Orf_Dict[chrom] = dict() 
	for pos in range(orf_start,orf_end+1):
		Pos_Orf_Dict[chrom][pos] = orf_name


################################################################################
Read_Starts_Dict = dict()
Read_Lengths_Dict = dict()

# now process the sam file
samFile = open(sys.argv[5],'r')
count = 0 
for line in samFile:
	if (line.startswith("@")):
		continue
	ts = line.strip().split()
	# check for unmapped reads
	if ((ts[2] == "*") or (ts[3] == "*")): 
		continue
	# check for unmapped reads
	if (int(ts[1]) & 0x0004):
		continue
	# check that mapping quality is > 10
	if (int(ts[4]) < 10):
		continue
    # check that the XM:i:<N>, if present, is <= 2
	mismatch = 0
	for i in range(11, len(ts)):
		if ("XM:i" in ts[i]):
			mismatch = int(re.compile("XM:i:([0-9]*)").search(ts[i]).group(1))
			break
	if (mismatch > 2):
		continue
	chrom = ts[2]
	if (chrom not in Pos_Orf_Dict): #ignore mitochondria, etc
		continue
	read_start = int(ts[3]) - 1 # into the 0-based coordinates
	# parse the cigar string to get the end		
	# find M for match
	ms = re.findall(r'([0-9]*)M',ts[5]) 
	ms = list(map(int, ms))
	Ms = sum(ms)
	#Ms = int(re.compile("([0-9]*)M").search(ts[5]).group(1))
	# find D for deletion 
	ds = re.findall(r'([0-9]*)D',ts[5]) 
	ds = list(map(int, ds))
	Ds = sum(ds)
	#if (re.compile("([0-9]*)D").search(ts[5]) != None):
	#	Ds = int(re.compile("([0-9]*)D").search(ts[5]).group(1))
	read_end = read_start + Ms + Ds -1
	read_len = Ms + Ds	
	# check that start is inside its gene orf
	if (read_start not in Pos_Orf_Dict[chrom]):
		continue
	gene = Pos_Orf_Dict[chrom][read_start]

	if (read_len == 28 or read_len == 25): # no need to shift reads, majority class is in frame 1
		pass
	elif (read_len == 29 or read_len == 31):  # majority class is in frame 3, shift forward by 1
		if (Orfs_Seq_Dict[gene]["strand"] == "-"):
			read_end -= 1
			read_len -= 1
		else:
			read_start += 1
			read_len -= 1
	elif (read_len == 32 or read_len == 27):  # majority class is in frame 2, shift back by 1
		if (Orfs_Seq_Dict[gene]["strand"] == "-"): 
			read_end += 1
			read_len += 1
		else:
			read_start -= 1
			read_len += 1
	else:
		continue

	
	# get the statistics for the read starts
	frame = -1
	if gene not in Read_Starts_Dict:
		Read_Starts_Dict[gene] = dict()
	if (Orfs_Seq_Dict[gene]["strand"] == "-"):
		frame = (Orfs_Seq_Dict[gene]["start"] - read_end)%3;
		if read_end not in Read_Starts_Dict[gene]:
			Read_Starts_Dict[gene][read_end] = 0
		Read_Starts_Dict[gene][read_end] = Read_Starts_Dict[gene][read_end] + 1
	else:	
		frame = (read_start - Orfs_Seq_Dict[gene]["start"])%3 
		if read_start not in Read_Starts_Dict[gene]:
			Read_Starts_Dict[gene][read_start] = 0
		Read_Starts_Dict[gene][read_start] = Read_Starts_Dict[gene][read_start] + 1
	if read_len not in Read_Lengths_Dict:
		Read_Lengths_Dict[read_len] = dict() 
		Read_Lengths_Dict[read_len][1] = 0	
		Read_Lengths_Dict[read_len][2] = 0	
		Read_Lengths_Dict[read_len][3] = 0	
	Read_Lengths_Dict[read_len][frame + 1] += 1
		#print("Stats for read")
		#print(ts)
		#print("read start: " + str(read_start) + " read length " + str(read_len))
		#print("gene " + gene + " gene_start " + str(Orfs_Seq_Dict[gene]["start"]) + " frame is " + str(frame))
		#print(Orfs_Seq_Dict[gene]["seq"])
		#print(Orfs_Seq_Dict[gene]["seq"][(read_start - Orfs_Seq_Dict[gene]["start"]):(read_start - Orfs_Seq_Dict[gene]["start"])+read_len])
	count +=1
	#if (count % 50000  == 0):
	#	print(count)
	#	break	
print(Read_Lengths_Dict)
for key in sorted(Read_Lengths_Dict.keys()):
	total = Read_Lengths_Dict[key][1] + Read_Lengths_Dict[key][2] + Read_Lengths_Dict[key][3]
	print("Length " + str(key) + " " + str(total) + " reads,  1: " + str(round(float(Read_Lengths_Dict[key][1])/total,2)) + " 2: " + str(round(float(Read_Lengths_Dict[key][2])/total,2)) + " 3: " + str(round(float(Read_Lengths_Dict[key][3])/total,2)))


outFile = open(sys.argv[6],'w')
for gene in Read_Starts_Dict:
	outFile.write(gene + "\n")
	outFile.write(Orfs_Seq_Dict[gene]["seq"] + "\n")
	for i in range(0,len(Orfs_Seq_Dict[gene]["seq"])):
		if (Orfs_Seq_Dict[gene]["strand"] == "-"):
			idx = Orfs_Seq_Dict[gene]["start"] - i
		else:
			idx = Orfs_Seq_Dict[gene]["start"] + i
		if idx in Read_Starts_Dict[gene]:
			outFile.write(str(Read_Starts_Dict[gene][idx]) + ",")
		else:
			outFile.write("0,")
	outFile.write("\n")
outFile.close()
