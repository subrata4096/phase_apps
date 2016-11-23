#!/usr/bin/python
import sys
import math
import os.path
#goldFile = "lulesh_dump_1_1_1.txt"
goldFile = ""


def calculateSpeedup(goldFile, fileToCompare):
	fp1 = open(goldFile, 'r')
	fp2 = open(fileToCompare, 'r')

	gold_instructions = 1
	gold_iterations = 1
	this_instructions = 1
	this_iterations = 1

	for line in fp1.readlines():
		line = line.strip()
		pos1 = line.find("Number of HW instructions Used")
		pos2 = line.find("Iteration count")
		if(pos1 != -1):
			fields = line.split("=")
			val = float(fields[1].strip())
			gold_instructions = val
			continue	
		if(pos2 != -1):
			fields = line.split("=")
			val = float(fields[1].strip())
			gold_iterations = val
			continue	

	for line in fp2.readlines():
		line = line.strip()
		pos1 = line.find("Number of HW instructions Used")
		pos2 = line.find("Iteration count")
		if(pos1 != -1):
			fields = line.split("=")
			val = float(fields[1].strip())
			this_instructions = val
			continue	
		if(pos2 != -1):
			fields = line.split("=")
			val = float(fields[1].strip())
			this_iterations = val
			continue	

	
	speedup = float(gold_instructions/this_instructions)
	return speedup,this_iterations

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
								goldFile = "lulesh_out_instruction_count_" + input_param_1 + "_" + input_param_2 + "_"  +  str(1) + "_" + str(-1) + "_" + str(1) + "_" + str(1) + "_" + str(1) + ".txt"
								fname = "lulesh_out_instruction_count_" + input_param_1 + "_" + input_param_2 + "_"  +  str(rep) + "_" + str(phase) + "_" + str(a1) + "_" + str(a2) + "_" + str(a3) + ".txt"
								if(False == os.path.isfile(fname)):
									continue
								speedup,iterations = calculateSpeedup(goldFile,fname)			
								#print str(phase) + "," + str(a1) + "," + str(a2) + "," + str(a3) + "," + str(avgErr)
								#print str(totalPhase) + "," + str(phase) + "," + str(a1) + "," + str(a2) + "," + str(a3) + "," + str(avgErrCount)
								print str(phase) + "," + input_param_1 + "," + input_param_2 + "," + str(a1) + "," + str(a2) + "," + str(a3) + "," + str(iterations) + "," + str(speedup)
