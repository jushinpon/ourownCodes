#!/bin/bash
#SBATCH -J CGDNA_example
#SBATCH -p Intel #Sumit Partition
#SBATCH -o test.out #Log file (maybe u don't need it)
#SBATCH -N 1#how many nodes u want 2 use
#SBATCH -w node[08]

export PATH=/opt/mpich_download/mpich-3.3.1/mpich_install/bin:$PATH
export LD_LIBRARY_PATH=/opt/mpich_download/mpich-3.3.1/mpich_install/lib:$LD_LIBRARY_PATH
/opt/lammps/lmp_mpi_202011261358 -sf omp -pk omp 32 -in Example.in
#mpirun -np 32 /opt/lammps/lmp_mpi_202009241029 -in Example.in  # The total cpu have to match with the nodes cpu quantity you had use

#squeue       # qstat
#scancel  1   # Delete JOB NO. 1 
#scontrol show job 1  # Check job NO.1 detail



################################################## DEBUG
#yum install dos2unix  #######   install this when you got as follows error messages
#sbatch: error: Batch script contains DOS line breaks (\r\n)
#sbatch: error: instead of expected UNIX line breaks (\n).
#dos2unix filename 
