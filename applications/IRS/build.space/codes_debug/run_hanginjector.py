#!/usr/gapps/asde/python/chaos_5_x86_64_ib/bin/python

import sys
import os

input_file = sys.argv[1]
f = open(input_file, 'r')
k = 0
for line in f:
	line = line.strip()
	fields = line.split(",")
	func_name = fields[0]
	call_num = fields[1]
        
	print "\n ***********************************************************"
	print "Case ", k, "Injecting for function: " , func_name 
	print "------------------------------------------------------------"

        os.system("rm -f /g/g90/mitra3/lock.irs")
	os.system("rm -f uniquefile*")

	syscall = "date; ./hang_ing.sh " + func_name + " " + call_num + " >> injector_old_125.old" 
	os.system(syscall)
	print " ***************************************************************\n"
        k = k + 1

f.close()
