#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=3950M
#SBATCH --job-name=EDIT
#SBATCH --account=def-rstull

ulimit -s unlimited
ulimit -v unlimited
ulimit -m unlimited

cd /home/tchui/scratch/CLIMATEX/PGW4ERA5
./process_era5.sh EDIT

