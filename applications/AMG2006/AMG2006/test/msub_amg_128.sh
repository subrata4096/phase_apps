#!/bin/csh
##### These lines are for Moab
#MSUB -A guests
#MSUB -l partition=vulcan
#MSUB -l nodes=2048
#MSUB -l walltime=00:15:00
#MSUB -q pbatch
#MSUB -V
#MSUB -m ae
#MSUB -o /p/lscratchv/mitra3/myjob-AVL-16384.out
#MSUB -e /p/lscratchv/mitra3/myjob-AVL-16384.out

##### These are shell commands
date
cd /g/g90/mitra3/work/case_studies/ROSS/rossnet-build/ross/models/phold
srun -n16384 ./phold --synch=3 --memory=4000
date
echo 'Done'
