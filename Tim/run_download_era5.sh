#!/bin/bash -l

if [ $# -lt 6 ]; then
    # TODO: print usage
    echo "Usage: run_ungrib SSYEAR SSMONTH SSDAY EEYEAR EEMONTH EEDAY RUNNAME"
    echo "RUNNAME IS OPTIONAL"
    exit 1
fi


SSYEAR=$1
SSMONTH=$2
SSDAY=$3
EEYEAR=$4
EEMONTH=$5
EEDAY=$6

RUNNAME=$7

user=`whoami`

ROOT_DIR=/scratch/${user}/CLIMATEX

if [ ! -z "$RUNNAME" ]; then
    ungrib_dir=${ROOT_DIR}/WPS_${SSYEAR}${SSMONTH}_${EEYEAR}${EEMONTH}_${user}_${RUNNAME}
else
    ungrib_dir=${ROOT_DIR}/WPS_${SSYEAR}${SSMONTH}_${EEYEAR}${EEMONTH}_${user}
fi

IBCS_DIR=${ungrib_dir}/ibcs/

echo ""
echo "==================="
echo ""
echo "Making directory ${ungrib_dir}"
echo ""
echo "==================="
echo ""

mkdir -p ${IBCS_DIR}
cd ${IBCS_DIR} || exit 1

smonths=( "9" "10" "11" "12" )
emonths=( "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" )
emonthsp1=( "1" )
eyearp1=$(( EEYEAR + 1 ))

download_year=$SSYEAR
for download_month in "${smonths[@]}"; do 
    if [[ ${download_month} == 9 ]]; then
        download_land="Y"
    else
        download_land="N"
    fi

    download_month=`printf "%02d\n" ${download_month}`


    mkdir ${download_year}_${download_month}
    cd ${download_year}_${download_month}

    download_directory=`pwd`
    echo $download_directory

    if [ -e IC.OK ]; then
        echo "IC.OK found in $download_directory --> doing nothing"
    else
        cp /home/${user}/climatex/submit_run_download_era5.sh .

        sed -i -e "s%DOWNLOAD_DIRECTORY%${download_directory}%g" submit_run_download_era5.sh
        sed -i -e "s%DOWNLOAD_YEAR%${download_year}%g" submit_run_download_era5.sh
        sed -i -e "s%DOWNLOAD_MONTH%${download_month}%g" submit_run_download_era5.sh
        sed -i -e "s%DOWNLOAD_LAND%${download_land}%g" submit_run_download_era5.sh
        sed -i -e "s%DJOB%D${download_year}${download_month}%g" submit_run_download_era5.sh
        
        sbatch submit_run_download_era5.sh
    
    fi
    
    
    cd ${IBCS_DIR}

done

download_land="N"
download_year=$EEYEAR
for download_month in ${emonths[@]}; do 
    download_month=`printf "%02d\n" ${download_month}`

    mkdir ${download_year}_${download_month}
    cd ${download_year}_${download_month}

    download_directory=`pwd`
    echo $download_directory


    if [ -e IC.OK ]; then
        echo "IC.OK found in $download_directory --> doing nothing"
    else
        cp /home/${user}/climatex/submit_run_download_era5.sh .

       sed -i -e "s%DOWNLOAD_DIRECTORY%${download_directory}%g" submit_run_download_era5.sh
        sed -i -e "s%DOWNLOAD_YEAR%${download_year}%g" submit_run_download_era5.sh
        sed -i -e "s%DOWNLOAD_MONTH%${download_month}%g" submit_run_download_era5.sh
        sed -i -e "s%DOWNLOAD_LAND%${download_land}%g" submit_run_download_era5.sh
        sed -i -e "s%DJOB%D${download_year}${download_month}%g" submit_run_download_era5.sh
        
        sbatch submit_run_download_era5.sh

    fi
    
   
    cd ${IBCS_DIR}

done


download_land="N"
download_year=$eyearp1
for download_month in ${emonthsp1[@]}; do 

    download_month=`printf "%02d\n" ${download_month}`

    mkdir ${download_year}_${download_month}
    cd ${download_year}_${download_month}

    download_directory=`pwd`
    echo $download_directory
    

    if [ -e IC.OK ]; then
        echo "IC.OK found in $download_directory --> doing nothing"
    else
        cp /home/${user}/climatex/submit_run_download_era5.sh .

        sed -i -e "s%DOWNLOAD_DIRECTORY%${download_directory}%g" submit_run_download_era5.sh
        sed -i -e "s%DOWNLOAD_YEAR%${download_year}%g" submit_run_download_era5.sh
        sed -i -e "s%DOWNLOAD_MONTH%${download_month}%g" submit_run_download_era5.sh
        sed -i -e "s%DOWNLOAD_LAND%${download_land}%g" submit_run_download_era5.sh
        sed -i -e "s%DJOB%D${download_year}${download_month}%g" submit_run_download_era5.sh
        
        sbatch submit_run_download_era5.sh
    
    fi
    
    cd ${IBCS_DIR}

done



