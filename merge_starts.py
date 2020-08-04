#!/usr/bin/python

#infiles = ["output/Hanson2018_bowtie_unique_25_32_starts", "output/Gerashchenko2014_bowtie_unique_25_32_starts", "output/Guydosh2014_bowtie_unique_25_32_starts", "output/Williams2014_bowtie_unique_25_32_starts", "output/Young2015_bowtie_unique_25_32_starts"]
#infiles = ["output/Hanson2018_bowtie_unique_25_32_starts", "output/Gerashchenko2014_bowtie_unique_25_32_starts", "output/Guydosh2014_bowtie_unique_25_32_starts", "output/Williams2014_bowtie_unique_25_32_starts", "output/Young2015_bowtie_unique_25_32_starts", "output/Wu2019_bowtie_unique_25_32_starts", "output/all_jan2014_bowtie_unique_25_to_32_starts"]
#outFile = open("output/combined_starts_no_wu_no_jan", "w")
#outFile = open("output/combined_starts_all_seven", "w")
infiles = ["output/all_jan2014_bowtie_28_starts", "output/Wu2019_bowtie_unique_28_29_starts", "output/Guydosh2014_bowtie_unique_28_29_starts", "output/Young2015_bowtie_unique_28_starts"]
outFile = open("output/combined_starts_four_sets_above_90", "w")
seq_dict = dict()
starts_dict = dict()
for f in infiles:
	inFile = open(f, "r")
	gene = inFile.readline().strip()
	while gene:
		seq = inFile.readline()
		starts_string = inFile.readline()
		current_gene = gene
		gene = inFile.readline().strip()

		ss = starts_string.strip()[:-1].split(",")
		if (current_gene not in seq_dict):
			seq_dict[current_gene] = seq
			starts_dict[current_gene] = ss
		else:
			new_starts = [int(starts_dict[current_gene][i]) + int(ss[i]) for i in range(len(ss))]
			starts_dict[current_gene] = new_starts
	inFile.close()
for gene in seq_dict:
	outString = gene + "\n"
	outString += seq_dict[gene] 
	outString += ",".join([str(x) for x in starts_dict[gene]]) + ",\n"
	outFile.write(outString)
outFile.close()
