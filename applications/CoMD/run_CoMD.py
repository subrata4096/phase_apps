#!/usr/bin/python
import sys
import os

total_phase_list = [4]
reps = [1,2,3,4,5]
input_unitcells = [10,11,12,15,16,17,20,21,22,25]
phase_list = [-1, 1, 2, 3, 4]
#phase_list = [1, 2, 3, 4, -1]
for rep in reps:
        for in_uc in input_unitcells:
		for phase in phase_list:
			for a1 in range(1,6):
				for a2 in range(1,6):
					for a3 in range(1,6):
						outfile = "comd_out_instruction_count_" + str(in_uc) + "_" + str(rep) + "_" + str(phase) + "_" + str(a1) + "_" + str(a2) + "_" + str(a3) + "_Sept14.txt" 
						command = "export APPROX_PHASE=" + str(phase) + ";export APPROX_PHASE_TOTAL=" + str(4) + ";export LA1=" + str(a1) + ";export LA2=" + str(a2) + ";export LA3=" + str(a3) + ";export UC=" + str(in_uc) + ";export REP=" + str(rep) + ";./CoMD/bin/CoMD-serial -x " + str(in_uc) + " -y " + str(in_uc) + " -z " + str(in_uc) + " >> " + outfile
						print command
						os.system(command)
