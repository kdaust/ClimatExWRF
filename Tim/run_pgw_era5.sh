#!/bin/bash -l
set -x

if [ $# -lt 6 ]; then
    # TODO: print usage
    echo "Usage: run_ungrib SSYEAR SSMONTH SSDAY EEYEAR EEMONTH EEDAY"
    exit 1
fi


SSYEAR=$1
SSMONTH=$2
SSDAY=$3
EEYEAR=$4
EEMONTH=$5
EEDAY=$6

RUNNAME=$7

skip=$8


user=`whoami`

endyearp1=$(( EEYEAR + 1 ))

ROOT_DIR=/scratch/${user}/CLIMATEX

if [ ! -z "$RUNNAME" ]; then
    ungrib_dir=${ROOT_DIR}/WPS_${SSYEAR}${SSMONTH}_${EEYEAR}${EEMONTH}_${user}_${RUNNAME}
    ungrib_dir_base=${ROOT_DIR}/WPS_${SSYEAR}${SSMONTH}_${EEYEAR}${EEMONTH}_${user}

else
    exit 1
fi


IBCS_DIR=${ungrib_dir_base}/ibcs

echo ""
echo "==================="
echo ""
echo "Making directory ${ungrib_dir}"
echo ""
echo "==================="
echo ""

mkdir ${ungrib_dir}

cd ${ungrib_dir} || exit 1

mkdir PGW4ERA5
cd PGW4ERA5

mkdir ${RUNNAME}

cd ${RUNNAME}

months=(${SSYEAR}_09 ${SSYEAR}_10 ${SSYEAR}_11 ${SSYEAR}_12 ${EEYEAR}_01 ${EEYEAR}_02 ${EEYEAR}_03 ${EEYEAR}_04 ${EEYEAR}_05 ${EEYEAR}_06 ${EEYEAR}_07 ${EEYEAR}_08 ${EEYEAR}_09 ${EEYEAR}_10 ${EEYEAR}_11 ${EEYEAR}_12)

code=0
for each_month in ${months[@]};
    
    do echo $each_month
    mkdir $each_month

    cd $each_month
    cp -v ~/climatex/PGW/* .
    cp -v ~/climatex/pgw_era5.sh .
    rm -fv ERA5*${each_month}*grib
    ln -s ${IBCS_DIR}/ERA5*${each_month}*grib .

    cp -v ~/climatex/submit_pgw_era5.sh ./submit_pgw_era5_${runname}${each_month}.sh
    sed -i s/PGWJOB/${RUNNAME}${each_month}/g ./submit_pgw_era5_${runname}${each_month}.sh
    sed -i s@PIDIR@"${ungrib_dir}/PGW4ERA5"@g ./submit_pgw_era5_${runname}${each_month}.sh
    yyyy=`echo ${each_month} | cut -d'_' -f1`
    mm=`echo ${each_month} | cut -d'_' -f2`
    sed -i s/YYYY/${yyyy}/g ./submit_pgw_era5_${runname}${each_month}.sh
    sed -i s/MM/${mm}/g ./submit_pgw_era5_${runname}${each_month}.sh
    sed -i s/SKIP/${skip}/g ./submit_pgw_era5_${runname}${each_month}.sh
   
    if [[ $RUNNAME == 2 ]]; then
        sed -i s/PGW/2/g ./pgw_era5.sh
    elif [[ $RUNNAME == 3 ]]; then
        sed -i s/PGW/3/g ./pgw_era5.sh
    else
        sed -i s/PGW/${RUNNAME}/g ./pgw_era5.sh
    fi

    if [[ -e  ../${RUNNAME}${each_month}.OK  ]]; then
        code=$(( code + 0 ))
    else
        sbatch ./submit_pgw_era5_${runname}${each_month}.sh
        code=$(( code + 1 ))
    fi
    cd ..

done

if [[ $code == 0 ]]; then
    touch PGW.OK 
fi







