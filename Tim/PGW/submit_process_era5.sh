#!/bin/bash
#SBATCH --time=11:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=3950M
#SBATCH --job-name=2021_10
#SBATCH --account=def-rstull

ulimit -s unlimited
ulimit -v unlimited
ulimit -m unlimited

cd /home/tchui/scratch/CLIMATEX/PGW4ERA5
./process_era5.sh 2021_10

