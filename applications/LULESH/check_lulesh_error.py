#!/usr/bin/python
import sys
import math
import os.path
#goldFile = "lulesh_dump_1_1_1.txt"
goldFile = "lulesh_dump_-1_4_1_1_1.txt"


def calculateError(goldFile, fileToCompare):
	fp1 = open(goldFile, 'r')
	fp2 = open(fileToCompare, 'r')

	goldValues = []
	approxValues = []

	for line in fp1.readlines():
		line = line.strip()
		fields = line.split("=")
		val = float(fields[1].strip())
		goldValues.append(val)

	for line in fp2.readlines():
		line = line.strip()
		fields = line.split("=")
		val = float(fields[1].strip())
		approxValues.append(val)

	totalDiff = 0.0
	numElems = len(goldValues)
	count = 0
	for i in range(numElems):
		goldVal = goldValues[i]
        	approxVal = approxValues[i]
        	if(math.fabs(goldVal) < 0.0001):
			goldVal = 0
		if(goldVal == 0):
			continue
		diff = math.fabs((goldValues[i] - approxValues[i])/float(goldValues[i]))
        	#print str(diff) + "   " + str(i) + "    " + str(goldValues[i]) + "    " + str(approxValues[i])
		totalDiff = totalDiff + diff
		count = count + 1

	avgDiff = totalDiff/float(numElems)
	avgDiffCount = totalDiff/float(count)

	#print "Approximation error : " + str(avgDiff) + " :  error based on count : " + str(avgDiffCount)
	return avgDiff,avgDiffCount

if __name__ == "__main__":
	#fileToCompare = sys.argv[1]
	#print fileToCompare
	#input_param_1 =  sys.argv[1].strip()
	#input_param_2 =  sys.argv[2].strip()
	#goldFile = "lulesh_dump_-1_4_1_1_1_" + input_param_1 + "_" + input_param_2 + "_1.txt"
	
	phase_list = [1, 2, 3, 4, -1]
	total_phase_list = [4]
	input_param_1_list = ["10","15","20","25","30","35","40"]
	input_param_2 = "8"
	for rep in range(1,6):
                for input_param_1 in input_param_1_list:
			for totalPhase in total_phase_list:
        			for phase in phase_list:
                			if(phase > totalPhase):
                        			continue
			
					for a1 in range(1,6):
						for a2 in range(1,6):
							for a3 in range(1,6):
								goldFile = "lulesh_dump_-1_4_1_1_1_" + input_param_1 + "_" + input_param_2 + "_1.txt"
								fname = "lulesh_dump_" + str(phase) + "_" + str(totalPhase) + "_" + str(a1) + "_" + str(a2) + "_" + str(a3) + "_" + input_param_1 + "_" + input_param_2 + "_" + str(rep) + ".txt"
								if(False == os.path.isfile(fname)):
									continue
								avgErr,avgErrCount = calculateError(goldFile,fname)			
								#print str(phase) + "," + str(a1) + "," + str(a2) + "," + str(a3) + "," + str(avgErr)
								#print str(totalPhase) + "," + str(phase) + "," + str(a1) + "," + str(a2) + "," + str(a3) + "," + str(avgErrCount)
								print str(phase) + "," + input_param_1 + "," + input_param_2 + "," + str(a1) + "," + str(a2) + "," + str(a3) + "," + str(avgErrCount)
