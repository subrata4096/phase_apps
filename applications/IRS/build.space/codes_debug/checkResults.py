#!/usr/gapps/asde/python/chaos_5_x86_64_ib/bin/python

import sys
import os
import test


def getFaultyMachine(theline):
	f = theline.split("=")
	return f[1]

def getfaultyTask(faultyMc, taskline):
	taskline = taskline.rstrip(",")
	D = dict(item.split(":") for item in taskline.split(","))
	faulty = D[faultyMc]
	#print faulty
	faulty_int = int(faulty)
	return faulty


def isTaskInRange(faultyTask,range):
	f = range.split("-")
	low = int(f[0])
	high = int(f[1])
	#print faultyTask
	task = int(faultyTask)
	if((low <= task) and (task <= high)):
		return True
	else:
		return False

def isTaskpresent(faultyTask,task):
        if(task.find("-") == -1):
                if(faultyTask == task):
                        return True
                else:
                        return False
        else:
                return isTaskInRange(faultyTask,task)

def parseLPTaskLine(faultyTask,theline):
	#print faultyTask, "   " , theline
	found = False
	f1s = theline.split("[")
	begPart = f1s[1]
	f2s = begPart.split("]")
	taskStrings = f2s[0]
	taskgroups = taskStrings.split(",")
	if(len(taskgroups) == 1):
		task = taskgroups[0]
		return isTaskpresent(faultyTask,task)
	else:
		for task in taskgroups:
			found = isTaskpresent(faultyTask,task)
			if(found):
				return True
	

	return False

				




input_file = sys.argv[1]
f = open(input_file, 'r')

#per case
faultyProcLine = ""
faultyTaskLineFound = False
proc_to_task_found = False
proc_to_task_linenum = -1
fTask = ""
identified_faulty_task = False
lp_taskline_found = False
lp_taskline_num = -1
listOfStateLineSeen = False
num_of_lp_states = 0
this_LP_correctly_detected = False

#global
line_num = 0
total_LP_task_count = 0
total_LP_correctly_detected = 0
total_false_positive = 0

for line in f:
	line_num += 1
	#print fTask
	line = line.strip()
	if(line.find("AutomaDeD started in MPI_Init") != -1):
		faultyProcLine = ""
		faultyTaskLineFound = False
		proc_to_task_found = False
		proc_to_task_linenum = -1
		fTask = ""
		identified_faulty_task = False
		lp_taskline_found = False
		lp_taskline_num = -1
		listOfStateLineSeen = False
		num_of_lp_states = 0
		this_LP_correctly_detected = False
		continue

	#if(line.find("Injecting in") != -1):
		#print line

	if(line.find("PIN_FAULTY_PROCESS=") != -1):
		faultyProcLine = line
		faultyTaskLineFound = True
		continue
        
	if(faultyProcLine):
		if(line.find("PROCS_TO_TASKS_MAPPING") != -1):
			if(faultyTaskLineFound):
				proc_to_task_found = True
				proc_to_task_linenum = line_num
				continue

		if(proc_to_task_found and (line_num == (proc_to_task_linenum + 1))):
			faultyMachine = getFaultyMachine(faultyProcLine)
			fTask = getfaultyTask(faultyMachine,line)
			identified_faulty_task = True
			continue
			
		if(line.find("Number of states with LP-tasks") != -1):
			lp_taskline_found = True
			lp_taskline_num = line_num
			total_LP_task_count += 1
			continue
		if((listOfStateLineSeen == False) and (line.find("STATE") != -1) and lp_taskline_found and (line_num > lp_taskline_num)):
			num_of_lp_states += 1
			lp = parseLPTaskLine(fTask,line)
			if(lp):
				this_LP_correctly_detected = True
				total_LP_correctly_detected += 1
				#print "MATCH"
		if(line.find("States:") != -1):
			listOfStateLineSeen = True
			if((this_LP_correctly_detected) and (num_of_lp_states > 1)):
				fp = float(float(num_of_lp_states-1)/float(num_of_lp_states))
				total_false_positive += fp
				#print "FP: " , total_false_positive
			faultyProcLine = False
			continue
	
f.close()

accuracy = float(total_LP_correctly_detected/float(total_LP_task_count))
falsePositive = float(total_false_positive/float(total_LP_task_count))
#falsePositive = float(total_false_positive/float(total_LP_correctly_detected))

print "Detection accuracy: " + str(accuracy)
print "False positive: " + str(falsePositive)

print "Total correctly detected: " + str(total_LP_correctly_detected) + " out of: " + str(total_LP_task_count)
print "Total False positive: " + str(int(total_false_positive))
