#!/bin/bash
#SBATCH --time=22:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=3950M
#SBATCH --constraint="[skylake|cascade]"
#SBATCH --job-name=UNGRIB
#SBATCH --account=def-rstull


ulimit -s unlimited
ulimit -v unlimited
ulimit -m unlimited

~/climatex/run_ungrib.sh $1 $2 $3 $4 $5 $6 $7



