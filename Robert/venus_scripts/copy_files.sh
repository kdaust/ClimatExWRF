#!/bin/bash
#PBS -o /home/rpayne/venus_scripts/logs/output_merge.log
#PBS -e /home/rpayne/venus_scripts/logs/errors_merge.log
#PBS -m ae -M gpayne1654@uvic.ca
#PBS -N copy_files

set -x  # Echo commands as they run
echo "Job started on $(hostname) at $(date)"

echo "Copying files"
cp /home/rpayne/data-rpayne/unproc/WRF-USASK/ctl/Q2/*.nc /home/rpayne/data-rpayne/unproc/WRF-USASK/ctl/Q2_merged/ > /home/rpayne/venus_scripts/logs/copy_files.log 2>&1

echo "Job ended at $(date)"