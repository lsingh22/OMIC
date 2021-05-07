#!/usr/bin/env bash
#SBATCH -p wacc                         
#SBATCH --job-name=OMIC              

#SBATCH --nodes=2                 
#SBATCH --time=0-00:10:00               
#SBATCH --cpus-per-task=1          

#SBATCH --output="gpu.out"             
#SBATCH --error="gpu.err"                                 
#SBATCH --gres=gpu:1

module load cuda
module load openmpi
module load netcdf

#SEG SCANS, DONE FOR BOTH 1 & 2 RANKS
#for i in {1..8}
#do
#      mpirun -np 2 ./OMIC gpu $((2**i)) 1 11
#done

#for j in {1..8}
#do
#      mpirun -np 2 ./OMIC gpu 256 $((2**j)) 11
#done

