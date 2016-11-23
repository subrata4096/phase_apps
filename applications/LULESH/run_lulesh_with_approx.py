#!/usr/bin/python
import subprocess
import os

outfilename_prefix = "./approx_logs/lulesh_approx_log_"
maxLvl = 7
for i in range(1,maxLvl):
	for j in range(1,maxLvl):
		for k in range(1,maxLvl):
                        outfilename = outfilename_prefix + str(i) + "_" + str(j) + "_" + str(k) + ".txt"
			command = "export LA1=" + str(i) + "; export LA2=" + str(j) + "; export LA3=" + str(k) + "; ./lulesh2.0 > " + outfilename 
			#print command
                        os.system(command)
                        iterHandle = subprocess.Popen(['grep','Iter',outfilename], stdout= subprocess.PIPE)
			itercount = iterHandle.communicate()[0]
			#iterHandle.wait() 
			toPrint = str(i) + "," + str(j) + "," + str(k) + "," + str(itercount)
			print toPrint
