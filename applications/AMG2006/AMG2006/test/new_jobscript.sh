#!/bin/sh
#MSUB -A asccasc
#MSUB -l partition=cab
#MSUB -l nodes=32
#MSUB -l walltime=00:10:00
#MSUB -q pbatch
#MSUB -o /g/g90/mitra3/work/LULESH/test_dump_512.txt
#MSUB -e /g/g90/mitra3/work/LULESH/test_dump_512.txt
 
export PIN_INJECT=1
export PIN_LOCK_FILE="/g/g90/mitra3/lock.lulesh$3"
export PIN_FAULTY_FUNC=$1
export PIN_FAULTY_EXEC=$2
export PIN_ENTRY_FUNC=myEntry
export PIN_EXIT_FUNC=myExit
export PIN_PROFILE_IMG=lulesh2.0
export AUT_USE_CALL_PATH=TRUE
export AUT_TIMEOUT=60
export LD_LIBRARY_PATH=/g/g90/mitra3/work/for_cab/automaded_LA_24Oct/install/lib:/usr/local/tools/fftw/lib:/usr/global/tools/adept/callpath/chaos_5_x86_64_ib/lib:/usr/global/tools/adept-utils/chaos_5_x86_64_ib/lib 
srun -n 512 /g/g90/mitra3/work/pin-2.12-58423-gcc.4.4.7-linux/pin -t /g/g90/mitra3/work/pin-2.12-58423-gcc.4.4.7-linux/source/tools/profile_MPI_frequency/obj-intel64/profile_MPI_frequency.so -- /g/g90/mitra3/work/LULESH/lulesh2.0 -s 5 -i 300 >> dump_512_$3.txt
