#!/bin/bash
#SBATCH --time=11:59:00
#SBATCH --ntasks-per-node=192
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=3950M
#SBATCH --job-name=JOBNAME
#SBATCH --account=rrg-rstull

ulimit -s unlimited
ulimit -v unlimited
ulimit -m unlimited

RUNDIR=`pwd`

./process_pcic.sh JOBNAME

