#!/bin/bash -x

export PIN_INJECT=1
export PIN_LOCK_FILE="/g/g90/mitra3/lock.amg"
#export PIN_FAULTY_FUNC="myFunction"
export PIN_FAULTY_FUNC=$1
export PIN_FAULTY_EXEC=$2

export PIN_ENTRY_FUNC=myEntry
export PIN_EXIT_FUNC=myExit
export PIN_PROFILE_IMG=amg2006
#export AUT_DUMP_MM=TRUE
export AUT_USE_CALL_PATH=TRUE
#export AUT_DO_NOT_DUMP=TRUE

export AUT_TIMEOUT=60

export PIN_TOOL=/g/g90/mitra3/work/pin-2.12-58423-gcc.4.4.7-linux/source/tools/profile_MPI_frequency/obj-intel64/profile_MPI_frequency.so

#export AUT_DO_NOT_EXIT=TRUE

#PIN=/g/g90/laguna/projects/PIN/pin-2.12-58423-gcc.4.4.7-linux/pin

#PIN_TOOL=/g/g90/laguna/projects/PIN/pin-2.12-58423-gcc.4.4.7-linux/source/tools/ManualExamples/obj-intel64/hang_injector.so
#PIN_TOOL=$PIN_HANGINJ_BASE/hang_injector.so

#MPI_PROG="/g/g90/mitra3/work/AMG2006/AMG2006/test/amg2006 -P 2 2 2"
MPI_PROG="/g/g90/mitra3/work/AMG2006/AMG2006/test/amg2006 -P 8 8 8"

#AUTOMADED=/g/g90/laguna/projects/TDPS_experiments/lib/libfsm.so
#AUTOMADED=/g/g90/mitra3/work/code/automaded_main/install/lib/libstracker.so
#AUTOMADED=/g/g90/mitra3/work/code/automaded/lib/libfsm.so

#LD_PRELOAD=$AUTOMADED srun -ppdebug -n 8 pin -t $PIN_TOOL -- $MPI_PROG
srun -ppdebug -O -N 16 -n 512 pin -t $PIN_TOOL -- $MPI_PROG
#srun -ppdebug -n 8 pin -t $PIN_TOOL -- $MPI_PROG
