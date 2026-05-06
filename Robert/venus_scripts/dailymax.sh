#!/bin/bash -l
#PBS -o /home/rpayne/venus_scripts/logs/output_dmax.log
#PBS -e /home/rpayne/venus_scripts/logs/errors_dmax.log
#PBS -m ae -M gpayne1654@uvic.ca
#PBS -N UBCWRF_dmax
#PBS -l select=1:ncpus=1

# Activate Conda properly
cd /home/rpayne/
source /home/rpayne/miniconda3/etc/profile.d/conda.sh
conda activate data-proc

echo "Calculating daily maximum WS10 for WRF"
cdo daymax -selname,WS10 \
    /home/rpayne/data-rpayne/unproc/WRF-UBC/WS10/wrfubc_ctl_WS10_198910_202409.nc \
    /home/rpayne/data-rpayne/unproc/WRF-UBC/WS10/wrfubc_ctl_WS10_daymax_198910_202409.nc \
    > /home/rpayne/venus_scripts/logs/dmax.log 2>&1

echo "Job ended at $(date)"
