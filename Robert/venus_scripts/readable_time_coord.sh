#!/bin/bash
#PBS -m ae -M gpayne1654@uvic.ca
#PBS -o /home/rpayne/venus_scripts/logs/output_time.log
#PBS -e /home/rpayne/venus_scripts/logs/errors_time.log
#PBS -N ubcwrf_tcoord
#PBS -l select=2:ncpus=1

set -x  # Echo commands as they run
echo "Job started on $(hostname) at $(date)"

# Activate my python environment
source /home/rpayne/miniconda3/etc/profile.d/conda.sh
conda activate data-proc

python /home/rpayne/venus_scripts/readable_time_coord.py > /home/rpayne/venus_scripts/logs/process_tcoord.log 2>&1

echo "Job ended at $(date)"
