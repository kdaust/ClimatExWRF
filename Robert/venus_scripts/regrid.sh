#!/bin/bash
#PBS -o /home/rpayne/venus_scripts/logs/output_ymean.log
#PBS -e /home/rpayne/venus_scripts/logs/errors_ymean.log
#PBS -m ae -M gpayne1654@uvic.ca
#PBS -N regrid_vas
#PBS -l nodes=venus01:ppn=1

# Activate Conda
source /home/rpayne/miniconda3/etc/profile.d/conda.sh
conda activate data-proc

echo "Regridding"
cdo remapbil,/home/rpayne/data-rpayne/unproc/WRF-UBC/griddes.txt \
    /home/rpayne/data-rpayne/unproc/ERA5-LAND/WS/ws_dailymean_199001_202412.nc \
    /home/rpayne/data-rpayne/unproc/ERA5-LAND/WS/ws_dailymean_regrid_199001_202412.nc \
    > /home/rpayne/venus_scripts/logs/regrid.log 2>&1

    