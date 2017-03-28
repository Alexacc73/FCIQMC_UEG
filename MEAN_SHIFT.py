fil = open("SHIFT_114SO_rs0.5_PGENwork.txt", "r")
print "FILENAME = ", fil.name

summ1 = 0.0
summ2 = 0.0
size = 0
shiftArray = []
projArray = []

for line in fil.readlines():
	size += 1
	shiftArray.append(float(line.split(" ")[0]  ))
	projArray.append(float(line.split(" ")[1]  ))

fil.close()

mean_sample = int(size*0.45)

for num in range(0, mean_sample ):
	summ1 += shiftArray[-num-1]
	summ2 += projArray[-num-1]

summ1 = summ1 / float(mean_sample)
summ2 = summ2 / float(mean_sample)

print "MEAN SHIFT       = ", summ1
print "MEAN PROJECTOR   = ", summ2
print "DIFFERENCE (P-S) = ",summ2 - summ1 


