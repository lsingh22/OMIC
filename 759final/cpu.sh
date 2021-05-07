#!/usr/bin/env bash
#SBATCH -p wacc                         
#SBATCH --job-name=OMIC              

#SBATCH --nodes=2                
#SBATCH --time=0-00:05:00               
#SBATCH --cpus-per-task=12          

#SBATCH --output="cpu.out"             
#SBATCH --error="cpu.err"                                 

module load openmpi
module load netcdf

#RANK SCAN
#for i in {1..16}
#do
#      mpirun -np $i --bind-to none ./OMIC cpu 256 256 11
#done

#SEGMENT SCAN, DONE AT 11 THREADS and 2 PROCS
#for i in {1..8}
#do
#      mpirun -np 2 --bind-to none ./OMIC cpu $((2**i)) 1 11
#done

#for j in {1..8}
#do
#      mpirun -np 2 --bind-to none ./OMIC cpu 256 $((2**j)) 11
#done

#THREAD SCAN
#for i in {1..32}
#do
#   mpirun -np 1 --bind-to none ./OMIC cpu 256 2 $i
#done

