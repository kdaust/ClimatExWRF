#!/bin/bash -l
set -x
if [ $# -lt 2 ]; then
    # TODO: print usage
    echo "Usage: step_02_prepoc_deltas.py GCM_NAME GWL"
    echo "OPTIONAL: RUNNAME"
    exit 1
fi

gcm_name=$1
GWL=$2
if [[ ${gcm_name} == "ACCESS-ESM1-5" ]]; then
    echo "Processing ${gcm_name}, GWL=${GWL}"

elif [[ ${gcm_name} == "BCC-CSM2-MR" ]]; then
    echo "Processing ${gcm_name}, GWL=${GWL}"

elif [[ ${gcm_name} == "CNRM-ESM2-1" ]]; then
    echo "Processing ${gcm_name}, GWL=${GWL}"

elif [[ ${gcm_name} == "CanESM5" ]]; then
    echo "Processing ${gcm_name}, GWL=${GWL}"

elif [[ ${gcm_name} == "EC-Earth3" ]]; then
    echo "Processing ${gcm_name}, GWL=${GWL}"

elif [[ ${gcm_name} == "GFDL-ESM4" ]]; then
    echo "Processing ${gcm_name}, GWL=${GWL}"

elif [[ ${gcm_name} == "GISS-E2-1-G" ]]; then
    echo "Processing ${gcm_name}, GWL=${GWL}"

elif [[ ${gcm_name} == "INM-CM5-0" ]]; then
    echo "Processing ${gcm_name}, GWL=${GWL}"

elif [[ ${gcm_name} == "IPSL-CM6A-LR" ]]; then
    echo "Processing ${gcm_name}, GWL=${GWL}"

elif [[ ${gcm_name} == "MIROC6" ]]; then
    echo "Processing ${gcm_name}, GWL=${GWL}"

elif [[ ${gcm_name} == "MPI-ESM1-2-HR" ]]; then
    echo "Processing ${gcm_name}, GWL=${GWL}"

elif [[ ${gcm_name} == "MRI-ESM2-0" ]]; then
    echo "Processing ${gcm_name}, GWL=${GWL}"

elif [[ ${gcm_name} == "UKESM1-0-LL" ]]; then
    echo "Processing ${gcm_name}, GWL=${GWL}"

else
    echo "not a valid GCM"
    exit 1
fi


if [[ ${GWL} != 2 ]] && [[ ${GWL} != 3 ]] && [[ ${GWL} != pi ]]; then
    echo "not a valid GWL"
    exit 1
fi



deltas_dir=DELTAS/${gcm_name}_${GWL}
regridded_dir=REGRIDDED/${gcm_name}_${GWL}
mkdir ${deltas_dir}

ln -s `pwd`/step_01_extract_deltas/Amon/${gcm_name}/*${GWL}_monthly_delta.nc ${deltas_dir}/
ln -s `pwd`/step_01_extract_deltas/Omon/${gcm_name}/*${GWL}_monthly_delta.nc ${deltas_dir}/
ln -s `pwd`/step_01_extract_deltas/SImon/${gcm_name}/*${GWL}_monthly_delta.nc ${deltas_dir}/

ln -s `pwd`/step_01_extract_deltas/Amon/${gcm_name}/*deg1_monthly.nc ${deltas_dir}/
ln -s `pwd`/step_01_extract_deltas/Omon/${gcm_name}/*deg1_monthly.nc ${deltas_dir}/
ln -s `pwd`/step_01_extract_deltas/SImon/${gcm_name}/*deg1_monthly.nc ${deltas_dir}/

era5_sample=`ls ERA5_heatdome_input | head -n1`
python3 step_02_preproc_deltas.py -g ${GWL} -i ${deltas_dir} -e ERA5_heatdome_input/${era5_sample} -o ${regridded_dir} regridding
