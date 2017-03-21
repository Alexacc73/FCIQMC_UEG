fil = open("SHIFT_66SO_rs0.5_BOOLFIX.txt", "r")
print "FILENAME = ", fil.name

summ1 = 0.0
summ2 = 0.0
summ3 = 0.0
size = 0
shiftArray = []
projNumerArray = []
WrefArray = []

for line in fil.readlines():
	size += 1
	shiftArray.append(float(line.split(" ")[0]  ))
	projNumerArray.append(float(line.split(" ")[2]  ))
	WrefArray.append(float(line.split(" ")[3]  ))

fil.close()

mean_sample = int(size*0.90)

for num in range(0, mean_sample ):
	summ1 += shiftArray[-num-1]
	summ2 += projNumArray[-num-1]
	summ3 += WrefArray[-num-1]

summ1 = summ1 / float(mean_sample)
summ2 = summ2 / float(mean_sample)
summ3 = summ3 / float(mean_sample)

print "MEAN SHIFT       = ", summ1
print "SEPARATED MEAN PROJECTOR   = ", (summ2/summ3)
print "DIFFERENCE (P-S) = ", (summ2/summ3) - summ1 


