#!/usr/bin/python
import sys

all_lines = dict()
genePDict = dict()

with_stops_File = open(sys.argv[1],'r')
for line in with_stops_File:
	ts = line.strip().split(",")
	ts.append("stops")
	all_lines[ts[0]] = ts
	genePDict[ts[0]] = float(ts[1])
with_stops_File.close()
with_stops_File = open(sys.argv[2],'r')
for line in with_stops_File:
	ts = line.strip().split(",")
	ts.append("speed")
	all_lines[ts[0]] = ts
	genePDict[ts[0]] = float(ts[1])

sorted_PVals = sorted(genePDict.items(), key=lambda x: x[1])

genesOutFile = open(sys.argv[3], "w")
for (gene, p_val) in sorted_PVals:
	genesOutFile.write(",".join(all_lines[gene]))
	genesOutFile.write("\n")
with_stops_File.close()
genesOutFile.close()

