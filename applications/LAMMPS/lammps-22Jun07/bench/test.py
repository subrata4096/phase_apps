#!/usr/gapps/asde/python/chaos_5_x86_64_ib/bin/python
import sys
import os

#fname=sys.argv[0]
#exe_count = sys.argv[1]

#function = fname.strip()
#count = exe_count.strip()

def runPin(f,c):
	function = f.strip()
	count = c.strip()

	com = "rm -f /g/g90/mitra3/lock.lammps;rm -f uniquefile*;srun -ppdebug -O -N 16 -n 512 pin -t /g/g90/mitra3/work/pin-2.12-58423-gcc.4.4.7-linux/source/tools/profile_MPI_frequency/obj-intel64/profile_MPI_frequency.so -- /g/g90/mitra3/work/LAMMPS/lammps-22Jun07/bench/lmp_g++ -var x 40 -var y 40 -var z 40 < in.lj"
	full="date;export PIN_INJECT=1;export PIN_LOCK_FILE=\"/g/g90/mitra3/lock.lammps\";export PIN_ENTRY_FUNC=myEntry;export PIN_EXIT_FUNC=myExit;export PIN_PROFILE_IMG=lmp_g++;export AUT_USE_CALL_PATH=TRUE;export AUT_TIMEOUT=60;"
	#full="date;export PIN_INJECT=1;export PIN_LOCK_FILE=\"/g/g90/mitra3/lock.lammps\";export PIN_ENTRY_FUNC=myEntry;export PIN_EXIT_FUNC=myExit;export PIN_PROFILE_IMG=lmp_g++;export AUT_TIMEOUT=60;"
	full = full + "export PIN_FAULTY_EXEC=" + count + ";export PIN_FAULTY_FUNC=" + function + ";" + com + " >> injector_old_512.old"
	print full
	os.system(full)
