#!/bin/sh
# PBS job script: runs em2 with MPI on one node. The number of MPI
# processes must equal nproc in input.f90.
#PBS -l nodes=1:ppn=4
ulimit -s unlimited
export OMP_NUM_THREADS=2
cd $PBS_O_WORKDIR
NPROCS=`wc -l <$PBS_NODEFILE`
mpirun -machinefile $PBS_NODEFILE -np $NPROCS ./em2
