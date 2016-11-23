#!/usr/bin/python
import sys
import os

#phase_list = [-1, 1, 2, 3, 4, 5, 6, 7, 8]
#phase_list = [4, 5, 6, 7, 8]
phase_list = [1, 2, 3, 4, -1]
#phase_list = [3]
#total_phase_list = [2, 4, 8]
total_phase_list = [4]
reps = [1,2,3,4,5]
input_r = [8,9,10,11,12,13,14]
input_s = [10,15,20,25,30,35,40]
for rep in reps:
	for in_r in input_r:
		for in_s in input_s:
			for totalPhase in total_phase_list:
				for phase in phase_list:
					if(phase > totalPhase):
						continue
					for a1 in range(1,6):
					#for a1 in range(4,5):
						for a2 in range(1,6):
							for a3 in range(1,6): 
								outfile = "lulesh_out_instruction_count_" + str(in_s) + "_" + str(in_r) + "_" + str(rep) + "_" + str(phase) + "_" + str(a1) + "_" + str(a2) + "_" + str(a3) + ".txt"
								command = "export APPROX_PHASE=" + str(phase) + ";export APPROX_PHASE_TOTAL=" + str(totalPhase) + ";export LA1=" + str(a1) + ";export LA2=" + str(a2) + ";export LA3=" + str(a3) + ";export S=" + str(in_s) + ";export R=" + str(in_r) + ";export REP=" + str(rep) + ";./lulesh2.0 -s " + str(in_s) + " -r " + str(in_r) + " >> " +  outfile
								print command
								os.system(command)
