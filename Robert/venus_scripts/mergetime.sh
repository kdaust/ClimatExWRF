#!/bin/bash
#PBS -m ae -M gpayne1654@uvic.ca
#PBS -o /home/rpayne/venus_scripts/logs/output_merge.log
#PBS -e /home/rpayne/venus_scripts/logs/errors_merge.log
#PBS -N wrf_merge
#PBS -l select=1:ncpus=1
#PBS -l walltime=250:00:00

# Activate Conda
source /home/rpayne/miniconda3/etc/profile.d/conda.sh
conda activate data-proc

cdo mergetime /home/rpayne/data-rpayne/unproc/WRF-UBC/T2/unmerged/*.nc \
    /home/rpayne/data-rpayne/unproc/WRF-UBC/T2/wrfubc_ctl_T2_198910_202409.nc \
    > /home/rpayne/venus_scripts/logs/mergetime.log 2>&1

# cdo mergetime /home/rpayne/data-rpayne/unproc/WRF-UBC/WS10/unmerged/*.nc \
#     /home/rpayne/data-rpayne/unproc/WRF-UBC/WS10/wrfubc_ctl_WS10_198910_202409.nc \
#     >> /home/rpayne/venus_scripts/logs/mergetime.log 2>&1

cdo mergetime /home/rpayne/data-rpayne/unproc/WRF-UBC/PSFC/unmerged/*.nc \
    /home/rpayne/data-rpayne/unproc/WRF-UBC/PSFC/wrfubc_ctl_PSFC_198910_202409.nc \
    >> /home/rpayne/venus_scripts/logs/mergetime.log 2>&1

cdo mergetime /home/rpayne/data-rpayne/unproc/WRF-UBC/Q2/unmerged/*.nc \
    /home/rpayne/data-rpayne/unproc/WRF-UBC/Q2/wrfubc_ctl_Q2_198910_202409.nc \
    >> /home/rpayne/venus_scripts/logs/mergetime.log 2>&1

cdo mergetime /home/rpayne/data-rpayne/unproc/WRF-UBC/RH/unmerged/*.nc \
    /home/rpayne/data-rpayne/unproc/WRF-UBC/RH/wrfubc_ctl_RH_198910_202409.nc \
    >> /home/rpayne/venus_scripts/logs/mergetime.log 2>&1

# cdo mergetime /home/rpayne/data-rpayne/unproc/WRF-UBC/PREC/unmerged/*.nc \
#     /home/rpayne/data-rpayne/unproc/WRF-UBC/PREC/wrfubc_ctl_PREC_198910_202409.nc \
#     >> /home/rpayne/venus_scripts/logs/mergetime.log 2>&1