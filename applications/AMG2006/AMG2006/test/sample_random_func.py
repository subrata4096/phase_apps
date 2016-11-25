#!/usr/gapps/asde/python/chaos_5_x86_64_ib/bin/python

import sys
import os
import random

input_file = sys.argv[1]
sample_count = int(sys.argv[2])
f = open(input_file, 'r')
lines = f.readlines()
name_dict = dict()
uniq_func_dict = dict()
for line in lines:
	line = line.strip()
	funcs_fields = line.split(":")
	info = funcs_fields[1]
	fields = info.split()
	print info
	funcs = fields[0].strip()
	#if((funcs.startswith("_") or (funcs.startswith(".")))):
	#	continue
	func_name = funcs
	call_num = int(fields[3].strip())
	uniq_func_dict[func_name] = call_num
f.close()

total_num_of_funcs = len(uniq_func_dict)
#print total_num_of_funcs

i = 0
for f in uniq_func_dict:
#        print f
	num = uniq_func_dict[f]
	if(num <= 1):
		continue
#	print i
	name_dict[i] = (f,num)
	i = i + 1
total_num_of_uniq_funcs = len(name_dict)
items = [x for x in range(total_num_of_uniq_funcs)]
#print i
#print items
#sys.exit(2)
random.shuffle(items)
selected_samples = random.sample(items,  sample_count)
for k in selected_samples:
	fname,callcount = name_dict[k]
	count = 1
	if(callcount < 100):
        	count = random.randint(2,callcount)
	else:
        	count = random.randint(2,100)

	op = fname + "," + str(count)
	print op
