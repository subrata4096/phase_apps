#!/usr/gapps/asde/python/chaos_5_x86_64_ib/bin/python

import sys
import os
import test

input_file = sys.argv[1]
f = open(input_file, 'r')
i = 1
for line in f:
	line = line.strip()
	fields = line.split(",")
	func_name = fields[0]
	call_num = fields[1]
        
	print "\n ***********************************************************"
	print "CASE:", i , " Injecting for function: " , func_name 
	print "------------------------------------------------------------"

        os.system("rm -f /g/g90/mitra3/lock.lammps")
	os.system("rm -f uniquefile*")

	#syscall = "date; ./hang_ing.sh " + func_name + " " + call_num + " >> injector_output_file.txt" 
	#syscall = "date; ./hang_ing.sh " + func_name + " " + call_num + " >> injector_output_file.txt" 
	#os.system(syscall)
	test.runPin(func_name,call_num)
	print " ***************************************************************\n"
       
        i = i + 1

f.close()
