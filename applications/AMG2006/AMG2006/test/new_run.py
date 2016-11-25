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
        
        s = 'export PIN_INJECT=1;export PIN_LOCK_FILE="/g/g90/mitra3/lock.amg' + str(k) + '";'
	s += "export PIN_FAULTY_FUNC=" + func_name + ";"
	s += 'export PIN_FAULTY_EXEC=' + call_num + ';'
	s += "export PIN_ENTRY_FUNC=myEntry;export PIN_EXIT_FUNC=myExit;"
	s += "export PIN_PROFILE_IMG=amg2006;"
	s += "export AUT_USE_CALL_PATH=TRUE;export AUT_TIMEOUT=60;"
	s += "export LD_LIBRARY_PATH=/usr/local/tools/mvapich2-gnu-1.7/lib:/g/g90/mitra3/work/for_cab/automaded_old/install/lib:/usr/local/tools/fftw/lib:/usr/global/tools/adept/callpath/chaos_5_x86_64_ib/lib:/usr/global/tools/adept-utils/chaos_5_x86_64_ib/lib;\n\n"
	s += "srun -O -N 8 -n 256 /g/g90/mitra3/work/pin-2.12-58423-gcc.4.4.7-linux/pin -t /g/g90/mitra3/work/pin-2.12-58423-gcc.4.4.7-linux/source/tools/profile_MPI_frequency/obj-intel64/profile_MPI_frequency.so -- /g/g90/mitra3/work/AMG2006/AMG2006/test/amg2006 -P 8 8 4 >> /g/g90/mitra3/work/AMG2006/AMG2006/test/injector_old_256.old;"
	s += 'date;'
	print s
	k += 1

f.close()
