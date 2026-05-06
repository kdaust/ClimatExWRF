#!/bin/bash
#PBS -m ae -M gpayne1654@uvic.ca
#PBS -o /home/rpayne/venus_scripts/logs/output_ubcwrf.log
#PBS -e /home/rpayne/venus_scripts/logs/errors_ubcwrf.log
#PBS -N process_UBCWRF
#PBS -l select=2:ncpus=1
#PBS -l walltime=250:00:00

set -x  # Echo commands as they run
echo "Job started on $(hostname) at $(date)"

# Activate my python environment
source /home/rpayne/miniconda3/etc/profile.d/conda.sh
conda activate data-proc

# Run the Python script and redirect output to the log file
python /home/rpayne/venus_scripts/process_WRF_UBC_data.py > /home/rpayne/venus_scripts/logs/process_ubcwrf.log 2>&1
echo "Job ended at $(date)"

# Run the Python script and redirect output to the log file
# python /home/rpayne/venus_scripts/process_WRF_UBC_PRECdata.py > /home/rpayne/venus_scripts/logs/process_ubcwrf_PREC.log 2>&1

echo "Job ended at $(date)"