#!/bin/bash
set -euo pipefail

# Resolve the scripts before changing to the requested data directory.
if [[ ${1:-} == "--help" || ${1:-} == "-h" ]]; then
    echo "Usage: bash process_pcic.sh [DATA_DIRECTORY]"
    echo "Defaults to the current directory. Downloads, outputs and markers go here."
    exit 0
fi
if (( $# > 1 )) || [[ ${1-} == "" && $# -gt 0 ]]; then
    echo "Usage: bash process_pcic.sh [DATA_DIRECTORY]" >&2
    exit 2
fi
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
data_dir=${1:-.}
mkdir -p -- "$data_dir"
cd -- "$data_dir"
data_dir=$(pwd -P)
echo "Data directory: $data_dir"

# Tune against available RAM and disk throughput; Python alone defaults to 1.
PCIC_WORKERS=${PCIC_WORKERS:-4}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
#
#

echo ""
echo ""
echo " ############# "
# Format: YYYY-MM-DD HH:MM:SS
current_time=$(date "+%Y-%m-%d %H:%M:%S")
echo $current_time
echo " ############ "
echo ""
echo ""

metgrid_month=metgrid_2021_11
metgrid_prevmonth=metgrid_2021_10


#### Hard part ---> handle WRFOUT files
if [[ -e WRFOUT.OK ]]; then
    echo "WRFOUT.OK found"
else
    
    # Previous month, containing hour 00 for current month
    /blue/data/WRF/kd_processed/rclone-v1.70.3-linux-amd64/rclone copy -vc climatex:Share_Data/WPS_202109_202212_howard/${metgrid_prevmonth}/WRFOUT --include=wrfout_d03_2021-11-01_00* ./
    
    # All other files
    /blue/data/WRF/kd_processed/rclone-v1.70.3-linux-amd64/rclone copy -vc climatex:Share_Data/WPS_202109_202212_howard/${metgrid_month}/WRFOUT --include=wrfout_d03* ./
    
    # NetCDF reads compressed inputs directly. Keep only one copy per hour
    # (native or _compressed); the processor rejects duplicates and gaps.
    # A fresh directory prevents stale hourly outputs entering this merge.
    pcic_workdir=$(mktemp -d ./pcic-hours.XXXXXX)
    python "$script_dir/output_pcic.py" d03 --workers "$PCIC_WORKERS" \
        --input-dir "$data_dir" --output-dir "$pcic_workdir"
    
    # Merge files
    cdo mergetime "$pcic_workdir"/wrfpcic* WRFPCIC_INTERIM.nc

    touch WRFOUT.OK
    echo ""
    echo ""
    echo " ############# "
    # Format: YYYY-MM-DD HH:MM:SS
    current_time=$(date "+%Y-%m-%d %H:%M:%S")
    echo $current_time
    echo " ############ "
    echo ""
    echo ""

fi

# Download already prepared monthly file for precip
# DOES NOT INCLUE HSNOWNC; ignore for now (or ask if HSOLIDNC is adequate)
/blue/data/WRF/kd_processed/rclone-v1.70.3-linux-amd64/rclone copy -vc climatex2:Share_Data/SUBSETTED/2022/metgrid_2021_11/ --include=COMPRESSED*RAIN* ./
nccopy -d0 -s COMPRESSED_RAIN_d03_metgrid_2021_11.nc rain_interim.nc


# Grab wrfuvic variables
if [[ -e WRFUVIC.OK ]]; then
    echo "WRFUVIC.OK found"
else
    
    # Hour 00
    # Actually don't include because then it will be 1 timestep longer than rain and pcic
    #rclone copy -vc climatex:Share_Data/WPS_202109_202212_howard/${metgrid_prevmonth}/WRFUVIC --include=wrfuvic_d03_2021-11-01_00* .
    
    # All other files
    /blue/data/WRF/kd_processed/rclone-v1.70.3-linux-amd64/rclone copy -vc climatex:Share_Data/WPS_202109_202212_howard/${metgrid_month}/WRFUVIC --include=wrfuvic_d03* ./
    
    # Grab only requested variables
    for wrfuvic in `ls wrfuvic*d03*`; do
        if [[ $wrfuvic ==  wrfuvic_d03_2021-11-01_00:00:00 ]]; then
            echo "Removing $wrfuvic because we want to start on hour 01"
            rm -fv $wrfuvic
            continue 
        fi
        cdo -selvar,Q2,T2,PSFC,U10,V10,GRDFLX,HFX,QFX,LH,SNOW,SNOWH,SNOWC  $wrfuvic smaller_$wrfuvic
        rm -fv $wrfuvic
    done
     
    # Merge files
    cdo mergetime smaller_* ./wrfuvic_interim.nc
   
    touch WRFUVIC.OK

    echo ""
    echo ""
    echo " ############# "
    # Format: YYYY-MM-DD HH:MM:SS
    current_time=$(date "+%Y-%m-%d %H:%M:%S")
    echo $current_time
    echo " ############ "
    echo ""

fi


cdo merge WRFPCIC_INTERIM.nc rain_interim.nc wrfuvic_interim.nc WRFPCIC_metgrid_2021_11.nc
nccopy -d1 -s WRFPCIC_metgrid_2021_11.nc COMPRESSED_WRFPCIC_metgrid_2021_11.nc
