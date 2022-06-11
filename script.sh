#!/bin/bash -l

#$ -cwd
#$ -q isis_4+_threads

n_proc=$7

#echo /usr/lib64/openmpi/bin/mpirun -np $n_proc bash -c "ulimit -s unlimited && ./main_async' $1 $2 $3 $4 $5 $6"


#this two are the main option you need|!!!!!!!!!!!!!!
#collect /usr/lib64/openmpi/bin/mpirun -np $n_proc -- bash -c "ulimit -s unlimited && ./main_async $1 $2 $3 $4 $5 $6"
mpirun --cpus-per-proc 1 --bind-to none -np $n_proc bash -c "ulimit -s unlimited && ./main_async $1 $2 $3 $4 $5 $6"







#/usr/lib64/openmpi/bin/mpirun -np $n_proc bash -c "ulimit -s unlimited && ./main_async $1 $2 $3 $4 $5 $6"