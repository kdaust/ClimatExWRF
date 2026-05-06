#!/bin/bash -l
#PBS -o /home/rpayne/venus_scripts/logs/output_dmean.log
#PBS -e /home/rpayne/venus_scripts/logs/errors_dmean.log
#PBS -m ae -M gpayne1654@uvic.ca
#PBS -N UBCWRF_dmean
#PBS -l select=2:ncpus=1

# Activate Conda properly
# source /home/rpayne/miniconda3/etc/profile.d/conda.sh
# conda activate data-proc

echo "Calculating daily mean for WRF"
# cdo daymean -selname,WS10 \
#     /home/rpayne/data-rpayne/unproc/WRF-UBC/WS10/wrfubc_ctl_WS10_198910_202409.nc \
#     /home/rpayne/data-rpayne/unproc/WRF-UBC/WS10/wrfubc_ctl_WS10_daymean_198910_202409.nc \
#     > /home/rpayne/venus_scripts/logs/dmean.log 2>&1

cdo daymean -selname,Q2 \
    /home/rpayne/data-rpayne/unproc/WRF-UBC/Q2/wrfubc_ctl_Q2_198910_202409.nc \
    /home/rpayne/data-rpayne/unproc/WRF-UBC/Q2/wrfubc_ctl_Q2_daymean_198910_202409.nc \
    >> /home/rpayne/venus_scripts/logs/dmean.log 2>&1

echo "Job ended at $(date)"