#!/bin/bash -x

export PIN_INJECT=1
export PIN_LOCK_FILE="/g/g90/mitra3/lock.irs"
#export PIN_FAULTY_FUNC="myFunction"
export PIN_FAULTY_FUNC=$1
export PIN_FAULTY_EXEC=$2

#export AUT_DO_NOT_EXIT=TRUE

#PIN=/g/g90/laguna/projects/PIN/pin-2.12-58423-gcc.4.4.7-linux/pin

#PIN_TOOL=/g/g90/laguna/projects/PIN/pin-2.12-58423-gcc.4.4.7-linux/source/tools/ManualExamples/obj-intel64/hang_injector.so
PIN_TOOL=$PIN_HANGINJ_BASE/hang_injector.so

MPI_PROG="/g/g90/mitra3/work/IRS/build.space/codes_opt/irs zrad3d -child_io_off -def NDOMS=8"

#AUTOMADED=/g/g90/laguna/projects/TDPS_experiments/lib/libfsm.so
#AUTOMADED=/g/g90/mitra3/work/code/automaded_main/install/lib/libstracker.so
#AUTOMADED=/g/g90/mitra3/work/code/automaded/lib/libfsm.so

#LD_PRELOAD=$AUTOMADED srun -ppdebug -n 8 pin -t $PIN_TOOL -- $MPI_PROG
srun -ppdebug -n 8 pin -t $PIN_TOOL -- $MPI_PROG
