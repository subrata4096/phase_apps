#!/usr/gapps/asde/python/chaos_5_x86_64_ib/bin/python

import sys
import os
import random
import subprocess

input_file = sys.argv[1]
sample_count = int(sys.argv[2])
f = open(input_file, 'r')
lines = f.readlines()
name_dict = dict()
task_mapping = dict()
uniq_func_dict = dict()
for line in lines:
	line = line.strip()
	funcs_fields = line.split(":")
        info  = ""
        if(len(funcs_fields) > 1):
                info = funcs_fields[1]
        else:
                info = line	
	fields = info.split(",")
	
	funcs = fields[0].strip()
	#if((funcs.startswith("_") or (funcs.startswith(".")))):
	#	continue
	func_name = funcs
	call_num = int(fields[1].strip())
	sysstr = "c++filt " + func_name
	#os.system(sysstr)
	proc = subprocess.Popen(sysstr, stdout=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()
	out = out.strip()
	if(out.startswith("_")):
		continue
	if(out.find("std::") != -1):
		continue
	if(out.find(".plt") != -1):
                continue
	uniq_func_dict[func_name] = call_num
	if(func_name in task_mapping):
		mapcount = task_mapping[func_name]
		task_mapping[func_name] = mapcount + 1
	else:
		task_mapping[func_name] = 1
f.close()

total_num_of_funcs = len(uniq_func_dict)
#print total_num_of_funcs

i = 0
for f in uniq_func_dict:
#        print f
        mapcount = task_mapping[f]
	if(mapcount < 50):
		continue
        num = uniq_func_dict[f]
        if(num <= 3):
                continue
#       print i
        name_dict[i] = (f,num)
        i = i + 1
total_num_of_uniq_funcs = len(name_dict)
items = [x for x in range(total_num_of_uniq_funcs)]
print i
popu_size=len(items)
#print items
#sys.exit(2)
random.shuffle(items)
selected_samples = random.sample(items,  sample_count)

for k in selected_samples:
        fname,callcount = name_dict[k]
        count = 1
        if(callcount < 20):
                count = random.randint(2,callcount)
        else:
                count = random.randint(2,20)

        op = fname + "," + str(count)
        print op

