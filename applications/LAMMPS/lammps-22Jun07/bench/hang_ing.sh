#!/bin/bash -x

export PIN_INJECT=1
export PIN_LOCK_FILE="/g/g90/mitra3/lock.lammps"
#export PIN_FAULTY_FUNC="myFunction"
export PIN_FAULTY_FUNC=$1
export PIN_FAULTY_EXEC=$2

#export AUT_DO_NOT_EXIT=TRUE

#PIN=/g/g90/laguna/projects/PIN/pin-2.12-58423-gcc.4.4.7-linux/pin

#PIN_TOOL=/g/g90/laguna/projects/PIN/pin-2.12-58423-gcc.4.4.7-linux/source/tools/ManualExamples/obj-intel64/hang_injector.so
PIN_TOOL=$PIN_HANGINJ_BASE/hang_injector.so

#MPI_PROG="/g/g90/mitra3/work/LAMMPS/lammps-22Jun07/bench/lmp_g++ -var x 40 -var y 40 -var z 80 < in.lj"  
MPI_PROG="/g/g90/mitra3/work/LAMMPS/lammps-22Jun07/bench/lmp_g++ -var x 40 -var y 40 -var z 80 << in.lj"  
#eval $MPI_PROG
#AUTOMADED=/g/g90/laguna/projects/TDPS_experiments/lib/libfsm.so
#AUTOMADED=/g/g90/mitra3/work/code/automaded_main/install/lib/libstracker.so
#AUTOMADED=/g/g90/mitra3/work/code/automaded/lib/libfsm.so

#LD_PRELOAD=$AUTOMADED srun -ppdebug -n 8 pin -t $PIN_TOOL -- $MPI_PROG
#y=eval $MPI_PROG
srun -ppdebug -n 16 pin -t $PIN_TOOL -- $MPI_PROG
#srun -ppdebug -n 16 pin -t /g/g90/mitra3/work/pin-2.12-58423-gcc.4.4.7-linux/source/tools/hang_injector/obj-intel64/hang_injector.so -- /g/g90/mitra3/work/LAMMPS/lammps-22Jun07/bench/lmp_g++ -var x 40 -var y 40 -var z 80 < in.lj

#python test.py
