#!/bin/bash -l
set -x 

if [ $# -lt 6 ]; then
    # TODO: print usage
    echo "Usage: run_ungrib SSYEAR SSMONTH SSDAY EEYEAR EEMONTH EEDAY"
    echo "OPTIONAL: RUNNAME"
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
slack=/home/${user}/climatex/slack.sh
username="$(whoami)@$(hostname -s)"
channel=climatex_runs
target='<!channel>'

MYEAR=$(( SSYEAR + 1 ))
PYEAR=$(( SSYEAR ))
AYEAR=$(( SSYEAR + 2 ))

ROOT_DIR=/scratch/${user}/CLIMATEX

if [ ! -z "$RUNNAME" ]; then
    wps_dir=${ROOT_DIR}/WPS_${SSYEAR}${SSMONTH}_${EEYEAR}${EEMONTH}_${user}_${RUNNAME}
    WPSDIR=WPS_${SSYEAR}${SSMONTH}_${EEYEAR}${EEMONTH}_${user}_${RUNNAME}
else
    wps_dir=${ROOT_DIR}/WPS_${SSYEAR}${SSMONTH}_${EEYEAR}${EEMONTH}_${user}
    WPSDIR=WPS_${SSYEAR}${SSMONTH}_${EEYEAR}${EEMONTH}_${user}
fi

source ~/WRFV4.3.3/source_wrf

cd $wps_dir || exit 1

list_of_metdir=`ls | grep metgrid_`

for each_dir in ${list_of_metdir[@]}; do 
    
    ssyear=`echo ${each_dir} | cut -d"_" -f2`
    eeyear=$ssyear
    ssmonth=`echo ${each_dir} | cut -d"_" -f3`

    if [ $ssmonth == 01 ]; then
        eemonth=02
    elif [ $ssmonth == 02 ]; then
        eemonth=03
    elif [ $ssmonth == 03 ]; then
        eemonth=04
    elif [ $ssmonth == 04 ]; then
        eemonth=05
    elif [ $ssmonth == 05 ]; then
        eemonth=06
    elif [ $ssmonth == 06 ]; then
        eemonth=07
    elif [ $ssmonth == 07 ]; then
        eemonth=08
    elif [ $ssmonth == 08 ]; then
        eemonth=09
    elif [ $ssmonth == 09 ]; then
        eemonth=10
    elif [ $ssmonth == 10 ]; then
        eemonth=11
    elif [ $ssmonth == 11 ]; then
        eemonth=12
    elif [ $ssmonth == 12 ]; then
        eemonth=01
        eeyear=$(( ssyear + 1 ))
    else
        echo "Invalid month"
        exit 1
    fi    
    
    metgrid_dir=${each_dir}


    # Don't do anything if month is not finished
    if [ ! -e ${wps_dir}/${metgrid_dir}/WRF/REAL.OK ]; then
        echo "${wps_dir}/${metgrid_dir}/WRF/REAL.OK missing; exiting"
        echo ""
        exit 0
    fi
    

        ### GZIP OF REMAP FILES ###
        squeue -u ${user} | grep ${user:0:2}Z${ssyear:2}${ssmonth}_${eeyear:2}${eemonth}
        jobrunning=$?
        
        # Don't do anything if a job is still running
        if [[ $jobrunning == 0 ]]; then
            echo "Job ${user:0:2}Z${ssyear:2}${ssmonth}_${eeyear:2}${eemonth} still running"
            echo ""  
        else
            if [ ! -e ${wps_dir}/${metgrid_dir}/WRF/RCLONE_GZIP.OK ]; then
                cp /home/${user}/climatex/submit_run_gzip_wrf.sh ${wps_dir}/${metgrid_dir}/WRF/
                sed -i -e "s/JOBNAME/${user:0:2}Z${ssyear:2}${ssmonth}_${eeyear:2}${eemonth}/g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_gzip_wrf.sh
                sed -i -e "s%RUNDIR%${wps_dir}/${metgrid_dir}/WRF%g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_gzip_wrf.sh
                sed -i -e "s/WPSDIR/${WPSDIR}/g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_gzip_wrf.sh
                sed -i -e "s/METGRIDDIR/${metgrid_dir}/g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_gzip_wrf.sh

                sbatch ${wps_dir}/${metgrid_dir}/WRF/submit_run_gzip_wrf.sh

            fi
        fi

done






