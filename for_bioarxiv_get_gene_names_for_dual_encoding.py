#!/usr/bin/python

import sys
import re
import random
# gene	stop	p_value	fraction_shifted	pattern_pvalue

random.seed(10)

###########################################################################
################################################################################
# create a map of ORF names
orfsFile = open(sys.argv[1],'r')
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

genesWithStartsFile = open(sys.argv[2],'r')
genesOutFile = open(sys.argv[3], "w")
for line in genesWithStartsFile:
	ts = line.strip().split(",")
	genesOutFile.write(ts[0] + "," + Orfs_Seq_Dict[ts[0]]["gene_name"] + "," + ts[1] + "," + ts[2] + "," + ts[3] + "," + str(round(float((ts[4])))) + "," + ts[5] + "\n")
genesOutFile.close()
