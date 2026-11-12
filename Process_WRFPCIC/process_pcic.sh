#!/bin/bash
set -euo pipefail
set -x
#conda init
#conda activate pgw
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
cp -v ../*py .
cp -v ../*csv .

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

############ Change util paths here ############

#rclone_exec=/blue/data/WRF/kd_processed/rclone-v1.70.3-linux-amd64/rclone
rclone_exec=/cvmfs/soft.computecanada.ca/gentoo/2023/x86-64-v3/usr/bin/rclone
rclone_args="--transfers 1 --multi-thread-streams 2"
nccopy_exec=/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcc12/netcdf/4.9.2/bin/nccopy
cdo_exec=/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/MPI/gcc12/openmpi4/cdo/2.2.2/bin/cdo
python_exec=~/miniconda_v2/envs/pgw/bin/python3 

############ Change month here (can do programmatically) #########

#metgrid_month=metgrid_2021_11
metgrid_month=$1

################################

# Find availability in inventory file
month_line=`grep $metgrid_month wrfout_inventory.csv`
prev_month_line=`grep -B 1 $metgrid_month wrfout_inventory.csv | head -n1`

# Current month
metgrid_month_server=`echo $month_line | cut -d':' -f1`
year_str=`echo $metgrid_month | cut -d'_' -f2`
month_str=`echo $metgrid_month | cut -d'_' -f3`
directory=`echo $month_line | cut -d'/' -f2`

echo $metgrid_month_server
echo $directory

# Previous month
metgrid_prevmonth_server=`echo $prev_month_line | cut -d':' -f1`
metgrid_prevmonth=`echo $prev_month_line | cut -d'/' -f3`

echo $metgrid_prevmonth_server
echo $metgrid_prevmonth

#### Hard part ---> handle WRFOUT files
if [[ -e WRFOUT.OK ]]; then
    echo "WRFOUT.OK found"
else
    
    # Previous month, containing hour 00 for current month
    ${rclone_exec} copy -vc ${metgrid_prevmonth_server}:Share_Data/${directory}/${metgrid_prevmonth}/WRFOUT --include=wrfout_d03_${year_str}-${month_str}-01_00* ${rclone_args} ./
    ls -lh  wrfout_d03_${year_str}-${month_str}-01_00* 
    # Try again if wrfout doesn't exist, use wrfrst instead
    if [[ $? != 0 ]]; then
        ${rclone_exec} copy -vc ${metgrid_prevmonth_server}:Share_Data/${directory}/${metgrid_prevmonth}/WRFOUT --include=wrfrst_d03_${year_str}-${month_str}-01_00* ${rclone_args} ./
    fi

    # All other files
    ${rclone_exec} copy -vc ${metgrid_month_server}:Share_Data/${directory}/${metgrid_month}/WRFOUT --include=wrfout_d03* ${rclone_args} ./
    
    # NetCDF reads compressed inputs directly. Keep only one copy per hour
    # (native or _compressed); the processor rejects duplicates and gaps.
    # A fresh directory prevents stale hourly outputs entering this merge.
    pcic_workdir=$(mktemp -d ./pcic-hours.XXXXXX)
    ${python_exec} "$script_dir/output_pcic.py" d03 --workers "$PCIC_WORKERS" \
        --input-dir "$data_dir" --output-dir "$pcic_workdir"
    
    # Merge files
    ${cdo_exec} mergetime "$pcic_workdir"/wrfpcic* WRFPCIC_INTERIM.nc

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
# DOES NOT INCLUE HSNOWNC; need to work with wrfhourly files

wy=`echo "${directory:11:4}"`
${rclone_exec} copy -vc climatex2:Share_Data/SUBSETTED/${wy}/${metgrid_month}/ --include=COMPRESSED*RAIN* ${rclone_args} ./
compressed_rain=`ls COMPRESSED*RAIN*`
$nccopy_exec -d0 -s $compressed_rain rain_interim.nc


# Grab wrfuvic variables
if [[ -e WRFUVIC.OK ]]; then
    echo "WRFUVIC.OK found"
else
    
    # Hour 00
    # Actually don't include because then it will be 1 timestep longer than rain and pcic
    #rclone copy -vc climatex:Share_Data/WPS_202109_202212_howard/${metgrid_prevmonth}/WRFUVIC --include=wrfuvic_d03_2021-11-01_00* .
    
    # All other files
    ${rclone_exec} copy -vc ${metgrid_month_server}:Share_Data/${directory}/${metgrid_month}/WRFUVIC --include=wrfuvic_d03* ${rclone_args} ./
    
    # Grab only requested variables
    for wrfuvic in `ls wrfuvic*d03*`; do
        if [[ $wrfuvic ==  wrfuvic_d03_${year_str}-${month_str}-01_00:00:00 ]]; then
            echo "Removing $wrfuvic because we want to start on hour 01"
            rm -fv $wrfuvic
            continue 
        fi
        ${cdo_exec} -selvar,Q2,T2,PSFC,U10,V10,GRDFLX,HFX,QFX,LH,SNOW,SNOWH,SNOWC  $wrfuvic smaller_$wrfuvic
        rm -fv $wrfuvic
    done
     
    # Merge files
    ${cdo_exec} mergetime smaller_* ./wrfuvic_interim.nc
   
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

# Adjust descriptions and units for hourly precip variables
${cdo_exec} -setattribute,HPRECIPNC@long_name="Hourly accumulated total precipitation (all phases grid scale + convective)",HPRECIPNC@units="mm",HSOLIDNC@long_name="Hourly accumulated solid precipitation (snow graupel hail)",HSOLIDNC@units="mm",HLIQUIDNC@long_name="Hourly accumulated liquid precipitation (rain)",HLIQUIDNC@units="mm" rain_interim.nc  rain_interim_edited.nc

# Download wrfhourly files to grab hourly snowfall (HSOLIDNC)
${rclone_exec} copy -vc ${metgrid_month_server}:Share_Data/${directory}/${metgrid_month}/WRFHOURLY  . --include=wrfhourly*d03* ${rclone_args}

# Extract hourly accumulated snowfall and merge into an interim file
for wrfhourly in `ls wrfhourly*`; do ${cdo_exec} -selvar,HSNOWNC -setattribute,HSNOWNC@long_name="Hourly accumulated snowfall",HSNOWNC@units="mm" $wrfhourly snow_${wrfhourly}; done; ${cdo_exec} mergetime snow_* snow_interim.nc

# Merge all files together and compress
${cdo_exec} merge WRFPCIC_INTERIM.nc rain_interim_edited.nc snow_interim.nc wrfuvic_interim.nc WRFPCIC_${metgrid_month}.nc
$nccopy_exec -d1 -s WRFPCIC_${metgrid_month}.nc COMPRESSED_WRFPCIC_${metgrid_month}.nc



