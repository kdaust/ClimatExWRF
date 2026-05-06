#!/bin/bash -l
#PBS -o /home/rpayne/venus_scripts/logs/output_monmean.log
#PBS -e /home/rpayne/venus_scripts/logs/errors_monmean.log
#PBS -m ae -M gpayne1654@uvic.ca
#PBS -N UBCWRF_monmean
#PBS -l select=2:ncpus=1

# Activate Conda properly
# source /home/rpayne/miniconda3/etc/profile.d/conda.sh
# conda activate data-proc

echo "Calculating monthly mean for WRF"

cdo ymonmean -selname,WS10 -seldate,2000-10-01,2015-09-30\
    /home/rpayne/data-rpayne/unproc/WRF-UBC/WS10/wrfubc_ctl_WS10_198910_202409.nc \
    /home/rpayne/data-rpayne/unproc/WRF-UBC/WS10/wrfubc_ctl_WS10_ymonmean_200010_201509.nc \
    > /home/rpayne/venus_scripts/logs/monmean.log 2>&1

cdo ymonmean -selname,WS10 -seldate,2000-10-01,2015-09-30\
    /home/rpayne/data-rpayne/unproc/WRF-UBC/WS10/wrfubc_ctl_WS10_daymax_198910_202409.nc \
    /home/rpayne/data-rpayne/unproc/WRF-UBC/WS10/wrfubc_ctl_WS10_daymax_ymonmean_200010_201509.nc \
    >> /home/rpayne/venus_scripts/logs/monmean.log 2>&1

cdo ymonmean -selname,Q2 -seldate,2000-10-01,2015-09-30\
    /home/rpayne/data-rpayne/unproc/WRF-UBC/Q2/wrfubc_ctl_Q2_198910_202409.nc \
    /home/rpayne/data-rpayne/unproc/WRF-UBC/Q2/wrfubc_ctl_Q2_ymonmean_200010_201509.nc \
    >> /home/rpayne/venus_scripts/logs/monmean.log 2>&1

echo "Job ended at $(date)"