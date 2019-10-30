#!/usr/bin/python
import sys

inFile = open(sys.argv[1],'r')
gene = inFile.readline()
num_windows = 0


f1_sum = 0 
f2_sum = 0 
f3_sum = 0 
total = 0

r1_total = 0
r1_f1_sum = 0
r1_f2_sum = 0
r1_f3_sum = 0


r2_total = 0
r2_200_sum = 0
r2_020_sum = 0
r2_002_sum = 0
r2_110_sum = 0
r2_101_sum = 0
r2_011_sum = 0

r3_total = 0
r3_300_sum = 0
r3_030_sum = 0
r3_003_sum = 0
r3_210_sum = 0
r3_201_sum = 0
r3_120_sum = 0
r3_021_sum = 0
r3_102_sum = 0
r3_012_sum = 0
r3_111_sum = 0

r4_total = 0
r4_400_sum = 0
r4_040_sum = 0
r4_004_sum = 0
r4_310_sum = 0
r4_301_sum = 0
r4_130_sum = 0
r4_031_sum = 0
r4_103_sum = 0
r4_013_sum = 0
r4_211_sum = 0
r4_121_sum = 0
r4_112_sum = 0
r4_220_sum = 0
r4_202_sum = 0
r4_022_sum = 0

while gene:
	seq = inFile.readline()
	starts_string = inFile.readline()
	ss = starts_string.split(",")
	gene_total_reads = 0
	
	for i in range(0, len(seq) - 3, 3): # i is in frame
		f1_prob = int(ss[i])
		f2_prob = int(ss[i+1])
		f3_prob = int(ss[i+2])
		all_f = f1_prob + f2_prob + f3_prob
		if (all_f ==0):
			continue
		total += 1
		f1_sum += f1_prob/all_f
		f2_sum += f2_prob/all_f
		f3_sum += f3_prob/all_f
		# collect individual statistics
		if (all_f == 1):
			r1_total += 1
			if (f1_prob == 1):
				r1_f1_sum += 1
			elif (f2_prob == 1):
				r1_f2_sum += 1
			else:
				r1_f3_sum += 1
		if (all_f == 2):
			r2_total += 1
			if (f1_prob ==2):
				r2_200_sum += 1
			elif (f2_prob == 2):
				r2_020_sum += 1
			elif (f3_prob == 2):
				r2_002_sum += 1
			elif (f1_prob == 1 and f2_prob == 1):
				r2_110_sum += 1
			elif (f1_prob == 1 and f3_prob == 1):
				r2_101_sum += 1
			else:
				r2_011_sum += 1
		if (all_f == 3):
			r3_total += 1
			if (f1_prob == 3):
				r3_300_sum += 1 
			elif (f2_prob == 3):
				r3_030_sum += 1 
			elif (f3_prob == 3):
				r3_003_sum += 1
			elif (f1_prob == 2 and f2_prob == 1):
				r3_210_sum += 1 
			elif (f1_prob == 2 and f3_prob == 1):
				r3_201_sum += 1 
			elif (f1_prob == 1 and f2_prob == 2):
				r3_120_sum += 1 
			elif (f2_prob == 2 and f3_prob == 1): 
				r3_021_sum += 1 
			elif (f1_prob == 1 and f3_prob == 2):
				r3_102_sum += 1 
			elif (f2_prob == 1 and f3_prob == 2):
				r3_012_sum += 1 
			elif (f1_prob == 1 and f2_prob == 1 and f3_prob == 1):
				r3_111_sum += 1 
		if (all_f == 4):
			r4_total += 1
			if (f1_prob == 4):
				r4_400_sum += 1 
			elif (f2_prob == 4):
				r4_040_sum += 1 
			elif (f3_prob == 4):
				r4_004_sum += 1 
			elif (f1_prob == 3 and f2_prob == 1):
				r4_310_sum += 1 
			elif (f1_prob == 3 and f3_prob == 1):
				r4_301_sum += 1 
			elif (f1_prob == 1 and f2_prob == 3):
				r4_130_sum += 1 
			elif (f2_prob == 3 and f3_prob == 1):
				r4_031_sum += 1 
			elif (f1_prob == 1 and f3_prob == 3):
				r4_103_sum += 1 
			elif (f2_prob == 1 and f3_prob == 3):
				r4_013_sum += 1 
			elif (f1_prob == 2 and f2_prob == 1):
				r4_211_sum += 1 
			elif (f1_prob == 1 and f2_prob == 2):
				r4_121_sum += 1 
			elif (f1_prob == 1 and f2_prob == 1):
				r4_112_sum += 1 
			elif (f1_prob == 2 and f2_prob == 2):
				r4_220_sum += 1 
			elif (f1_prob == 2 and f2_prob == 0):
				r4_202_sum += 1 
			elif (f1_prob == 0 and f2_prob == 2):
				r4_022_sum += 1 

	gene = inFile.readline()
inFile.close()
f1_prob = f1_sum/total
f2_prob = f2_sum/total
f3_prob  = f3_sum/total
print("F1 average probability is " + str(round(f1_prob, 3)))
print("F2 average probability is " + str(round(f2_prob, 3)))
print("F3 average probability is " + str(round(f3_prob, 3)))

total_obs = 0
total_pred = 0
print("\nR1\tObs\tPred")
obs = r1_f1_sum/r1_total
pred = f1_prob
total_obs += obs
total_pred += pred 
print("100" + "\t" + str(round(obs,3)) + "\t" + str(round(pred, 3)))
obs = r1_f2_sum/r1_total
pred = f2_prob
total_obs += obs
total_pred += pred 
print("010" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r1_f3_sum/r1_total
pred = f3_prob
total_obs += obs
total_pred += pred 
print("001" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
print("total" + "\t" + str(round(total_obs, 3)) + "\t" + str(round(total_pred, 3)))

total_obs = 0
total_pred = 0
print("\nR2\tObs\tPred")
obs = r2_200_sum/r2_total
pred = f1_prob ** 2 
total_obs += obs
total_pred += pred 
print("200" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r2_020_sum/r2_total
pred = f2_prob ** 2 
total_obs += obs
total_pred += pred 
print("020" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r2_002_sum/r2_total
pred = f3_prob ** 2
total_obs += obs
total_pred += pred 
print("002" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r2_110_sum/r2_total
pred = 2.0 * f1_prob * f2_prob
total_obs += obs
total_pred += pred 
print("110" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r2_101_sum/r2_total
pred = 2.0 * f1_prob * f3_prob 
total_obs += obs
total_pred += pred 
print("101" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r2_011_sum/r2_total
pred = 2.0 * f2_prob * f3_prob
total_obs += obs
total_pred += pred 
print("011" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
print("total" + "\t" + str(round(total_obs, 3)) + "\t" + str(round(total_pred, 3)))


total_obs = 0
total_pred = 0
print("\nR3\tObs\tPred")
obs = r3_300_sum/r3_total
pred = f1_prob ** 3
total_obs += obs
total_pred += pred 
print("300" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r3_030_sum/r3_total
pred = f2_prob ** 3
total_obs += obs
total_pred += pred 
print("030" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r3_003_sum/r3_total
pred = f3_prob ** 3 
total_obs += obs
total_pred += pred 
print("003" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r3_210_sum/r3_total
pred = (f1_prob ** 2) * f2_prob * 3
total_obs += obs
total_pred += pred 
print("210" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r3_201_sum/r3_total
pred = (f1_prob ** 2) * f3_prob * 3
total_obs += obs
total_pred += pred 
print("201" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r3_120_sum/r3_total 
pred = f1_prob * (f2_prob ** 2) * 3 
total_obs += obs
total_pred += pred 
print("120" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r3_021_sum/r3_total
pred = (f2_prob ** 2) * f3_prob * 3 
total_obs += obs
total_pred += pred 
print("021" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r3_102_sum/r3_total 
pred = f1_prob * (f3_prob ** 2) * 3 
total_obs += obs
total_pred += pred 
print("102" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r3_012_sum/r3_total
pred = f2_prob * (f3_prob ** 2) * 3 
total_obs += obs
total_pred += pred 
print("012" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r3_111_sum/r3_total
pred = f1_prob * f2_prob * f3_prob * 6 
total_obs += obs
total_pred += pred 
print("111" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
print("total" + "\t" + str(round(total_obs, 3)) + "\t" + str(round(total_pred, 3)))

exit(0)
total_obs = 0
total_pred = 0
f1_prob = 0.33333333
f2_prob = 0.33333333
f3_prob = 0.33333333
print("\nR4\tObs\tPred")
obs = r4_400_sum/r4_total 
pred = f1_prob ** 4
total_obs += obs
total_pred += pred 
print("400" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r4_040_sum/r4_total 
pred = f2_prob ** 4
total_obs += obs
total_pred += pred 
print("040" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r4_004_sum/r4_total 
pred = f3_prob ** 4
total_obs += obs
total_pred += pred 
print("004" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r4_310_sum/r4_total 
pred = (f1_prob ** 3) * f2_prob * 4
total_obs += obs
total_pred += pred 
print("310" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r4_301_sum/r4_total 
pred = (f1_prob ** 3) * f3_prob * 4
total_obs += obs
total_pred += pred 
print("301" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r4_130_sum/r4_total 
pred = f1_prob * (f2_prob ** 3) * 4
total_obs += obs
total_pred += pred 
print("130" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r4_031_sum/r4_total 
pred = (f2_prob ** 3) * f3_prob * 4
total_obs += obs
total_pred += pred 
print("031" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r4_103_sum/r4_total 
pred = f1_prob * (f3_prob ** 3) * 4
total_obs += obs
total_pred += pred 
print("103" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r4_013_sum/r4_total 
pred = f2_prob * (f3_prob ** 3) * 4
total_obs += obs
total_pred += pred 
print("013" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r4_211_sum/r4_total 
pred = (f1_prob ** 2) * f2_prob * f3_prob * 12
total_obs += obs
total_pred += pred 
print("211" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r4_121_sum/r4_total 
pred = f1_prob * (f2_prob ** 2) * f3_prob * 12
total_obs += obs
total_pred += pred 
print("121" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r4_112_sum/r4_total 
pred = f1_prob * f2_prob * (f3_prob ** 2) * 12
total_obs += obs
total_pred += pred 
print("112" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r4_220_sum/r4_total 
pred = (f1_prob ** 2) * (f2_prob ** 2) * 12
total_obs += obs
total_pred += pred 
print("220" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r4_202_sum/r4_total 
pred = (f1_prob ** 2) * (f3_prob ** 2) * 12
total_obs += obs
total_pred += pred 
print("202" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
obs = r4_022_sum/r4_total 
pred = (f2_prob ** 2) * (f3_prob ** 2) * 12
total_obs += obs
total_pred += pred 
print("022" + "\t" + str(round(obs, 3)) + "\t" + str(round(pred, 3)))
print("total" + "\t" + str(round(total_obs, 3)) + "\t" + str(round(total_pred, 3)))
