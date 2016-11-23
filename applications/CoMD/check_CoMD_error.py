#!/usr/bin/python
import sys
import math

#goldFile = "output_CoMD_01_01_01.txt"
#goldFile = "error_profile_comd/output_CoMD_-1_04_01_01_01.txt"
#goldFile = "output_CoMD_-1_04_01_01_01.txt"
goldFile = "output_CoMD_01_04_01_01_01.txt"


def calculateError(goldFile, fileToCompare):
        fp1 = open(goldFile, 'r')
        fp2 = open(fileToCompare, 'r')

	goldValues = []
	approxValues = []

	for line in fp1.readlines():
		line = line.strip()
		fields = line.split(",")
		val1 = float(fields[1].strip())
		val2 = float(fields[2].strip())
		val3 = float(fields[3].strip())
		val4 = float(fields[4].strip())
		val5 = float(fields[5].strip())
		vals = [val1,val2,val3,val4,val5]
		goldValues.append(vals)

	for line in fp2.readlines():
		line = line.strip()
		fields = line.split(",")
		val1 = float(fields[1].strip())
		val2 = float(fields[2].strip())
		val3 = float(fields[3].strip())
		val4 = float(fields[4].strip())
		val5 = float(fields[5].strip())
		vals = [val1,val2,val3,val4,val5]
		approxValues.append(vals)

	totalDiff = 0.0
	numElems = len(goldValues)
	count = 0
	for i in range(numElems):
		goldVals = goldValues[i]
        	approxVals = approxValues[i]
		for k in range(5):
			goldVal = goldVals[k]
			approxVal = approxVals[k]
        		if(math.fabs(goldVal) < 0.00000001):
				goldVal = 0
			if(goldVal == 0):
				continue
			diff = math.fabs((goldVal - approxVal)/float(goldVal))
        		#print str(diff) + "   " + str(i) + "    " + str(goldValues[i]) + "    " + str(approxValues[i])
			totalDiff = totalDiff + diff
			count = count + 1

	avgDiff = totalDiff/float(numElems*5)
	avgDiffCount = totalDiff/float(count)

	#print "Approximation error : " + str(avgDiff) + " :  error based on count : " + str(avgDiffCount)
	return avgDiff, avgDiffCount

if __name__ == "__main__":
	#fileToCompare = sys.argv[1]
	#print fileToCompare
        for phase in [1,2,3,4,-1]:
                for a1 in range(1,6):
                        for a2 in range(1,6):
                                for a3 in range(1,6):
                                        #fname = "output_CoMD_" + str(phase) + "_4_" + str(a1) + "_" + str(a2) + "_" + str(a3) + ".txt"
                                        #fname = 'error_profile_comd/output_CoMD_%02d_04_%02d_%02d_%02d.txt' % (phase, a1, a2, a3)
                                        fname = 'output_CoMD_%02d_04_%02d_%02d_%02d.txt' % (phase, a1, a2, a3)
                                        avgErr,avgErrCount = calculateError(goldFile,fname)
                                        #print str(phase) + "," + str(a1) + "," + str(a2) + "," + str(a3) + "," + str(avgErr)
                                        print str(phase) + "," + str(a1) + "," + str(a2) + "," + str(a3) + "," + str(avgErrCount)
