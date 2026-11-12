#!/bin/bash  -l

isleap() {
   year=$1
   (( !(year % 4) && ( year % 100 || !(year % 400) ) )) &&
      echo "leap" || echo "nleap"
}


if [ $# -lt 2 ]; then
    # TODO: print usage
    echo "Usage: step_03_apply_to_era_pres.sh GCM_NAME GWL"
    echo "OPTIONAL: RUNNAME"
    exit 1
fi

gcm_name=$1
GWL=$2
yyyy=$3
mm=$4

source ~/miniconda_v2/bin/activate
export PATH=~/miniconda_v2/bin:$PATH
conda activate pgw

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

elif [[ ${gcm_name} == "ENS" ]]; then
    echo "Processing ${gcm_name}, GWL=${GWL}"

else
    echo "not a valid GCM"
    exit 1
fi


if [[ ${GWL} != 2 ]] && [[ ${GWL} != 3 ]] && [[ ${GWL} != pi ]]; then
    echo "not a valid GWL"
    exit 1
fi



regridded_dir=~/REGRIDDED/${gcm_name}_${GWL}
final_dir=../FINAL
mkdir ${final_dir}

if [[ $mm == 01 ]] ||  [[ $mm == 03 ]] || [[ $mm == 05 ]] || [[ $mm == 07 ]] || [[ $mm == 08 ]] || [[ $mm == 10 ]] || [[ $mm == 12 ]]; then
    dd=31
elif [[ $mm == 02 ]]; then
    res=`isleap $yyyy`
    if [[ ${res} == "leap" ]]; then
        dd=29
    else
        dd=28
    fi
else
    dd=30
fi


python3 step_03_apply_to_era_pres_v2.py -i ${yyyy}_${mm} -g ${GWL} -d ${regridded_dir} -o ${final_dir} -f ${yyyy}${mm}0100 -l ${yyyy}${mm}${dd}21


