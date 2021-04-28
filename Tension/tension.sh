#!/bin/bash
#SBATCH -J grain_tension
#SBATCH -p against_time #Sumit Partition
#SBATCH -o grain_tension.out #Log file (maybe u don't need it)
#SBATCH -N 1 #how many nodes u want 2 use
#SBATCH -w node[03]


#export PATH=/opt/mpich_download/mpich-3.3.1/mpich_install/bin:$PATH
#export LD_LIBRARY_PATH=/opt/mpich_download/mpich-3.3.1/mpich_install/lib:$LD_LIBRARY_PATH
export PATH=/opt/MVAPICH2/mvapich2-2.3.4/mpich_install/bin:$PATH
export LD_LIBRARY_PATH=/opt/MVAPICH2/mvapich2-2.3.4/mpich_install/lib:$LD_LIBRARY_PATH

mpirun -np 16 /opt/lammps/lmp_mpi_202009241712 -in Tension.in 

#squeue       # qstat
#scancel  1   # Delete JOB NO. 1 
#scontrol show job 1  # Check job NO.1 detail