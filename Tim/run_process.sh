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

ROOT_DIR=~/project/CLIMATEX

if [ ! -z "$RUNNAME" ]; then
    wps_dir=${ROOT_DIR}/WPS_${SSYEAR}${SSMONTH}_${EEYEAR}${EEMONTH}_${user}_${RUNNAME}
    WPSDIR=WPS_${SSYEAR}${SSMONTH}_${EEYEAR}${EEMONTH}_${user}_${RUNNAME}
else
    wps_dir=${ROOT_DIR}/WPS_${SSYEAR}${SSMONTH}_${EEYEAR}${EEMONTH}_${user}
    WPSDIR=WPS_${SSYEAR}${SSMONTH}_${EEYEAR}${EEMONTH}_${user}
fi

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
    if [ ! -e ${wps_dir}/${metgrid_dir}/WRF/WRF.OK ]; then
        echo "${wps_dir}/${metgrid_dir}/WRF/WRF.OK missing; exiting"
        echo ""
        exit 0
    fi

    # If all processing and transferring is finished, move on
    if [ ! -e ${wps_dir}/${metgrid_dir}/WRF/RCLONE.OK ]; then
        
        ### COMPRESSION OF WRFOUT FILES ###
        squeue -u ${user} | grep ${user:0:2}C${ssyear:2}${ssmonth}_${eeyear:2}${eemonth}
        jobrunning=$?
        
        # Don't do anything if a job is still running
        if [[ $jobrunning == 0 ]]; then
            echo "Job ${user:0:2}C${ssyear:2}${ssmonth}_${eeyear:2}${eemonth} still running"
            echo ""  
        else
            if [ ! -e ${wps_dir}/${metgrid_dir}/WRF/RCLONE_COMPRESS.OK ]; then
                cp /home/${user}/climatex/submit_run_compress_wrf.sh ${wps_dir}/${metgrid_dir}/WRF/
                cp /home/${user}/climatex/run_compress_wrf.sh ${wps_dir}/${metgrid_dir}/WRF/
                sed -i -e "s/JOBNAME/${user:0:2}C${ssyear:2}${ssmonth}_${eeyear:2}${eemonth}/g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_compress_wrf.sh
                sed -i -e "s%RUNDIR%${wps_dir}/${metgrid_dir}/WRF%g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_compress_wrf.sh
                sed -i -e "s/WPSDIR/${WPSDIR}/g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_compress_wrf.sh
                sed -i -e "s/METGRIDDIR/${metgrid_dir}/g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_compress_wrf.sh

                sbatch ${wps_dir}/${metgrid_dir}/WRF/submit_run_compress_wrf.sh

            fi
        fi

        ### REMAPPING OF WRFUVIC FILES ###
        squeue -u ${user} | grep ${user:0:2}V${ssyear:2}${ssmonth}_${eeyear:2}${eemonth}
        jobrunning=$?
        
        if [[ $jobrunning == 0 ]]; then
            echo "Job ${user:0:2}V${ssyear:2}${ssmonth}_${eeyear:2}${eemonth} still running"
            echo ""
        else
            if [ ! -e ${wps_dir}/${metgrid_dir}/WRF/RCLONE_WRFUVIC.OK ]; then
                cp /home/${user}/climatex/remapping/submit_run_regrid_wrfuvic.sh ${wps_dir}/${metgrid_dir}/WRF/
                cp /home/${user}/climatex/remapping/run_regrid_wrfuvic.sh ${wps_dir}/${metgrid_dir}/WRF/
                sed -i -e "s/JOBNAME/${user:0:2}V${ssyear:2}${ssmonth}_${eeyear:2}${eemonth}/g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_regrid_wrfuvic.sh
                sed -i -e "s%RUNDIR%${wps_dir}/${metgrid_dir}/WRF%g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_regrid_wrfuvic.sh
                sed -i -e "s/WPSDIR/${WPSDIR}/g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_regrid_wrfuvic.sh
                sed -i -e "s/METGRIDDIR/${metgrid_dir}/g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_regrid_wrfuvic.sh

                sbatch ${wps_dir}/${metgrid_dir}/WRF/submit_run_regrid_wrfuvic.sh
            fi
        fi

        ### REMAPPING OF WRFHOURLY FILES ###
        squeue -u ${user} | grep ${user:0:2}H${ssyear:2}${ssmonth}_${eeyear:2}${eemonth}
        jobrunning=$?
        
        if [[ $jobrunning == 0 ]]; then 
            echo "Job ${user:0:2}H${ssyear:2}${ssmonth}_${eeyear:2}${eemonth} still running"
            echo ""  
        else
            if [ ! -e ${wps_dir}/${metgrid_dir}/WRF/RCLONE_WRFHOURLY.OK ]; then
                cp /home/${user}/climatex/remapping/submit_run_regrid_wrfhourly.sh ${wps_dir}/${metgrid_dir}/WRF/
                cp /home/${user}/climatex/remapping/run_regrid_wrfhourly.sh ${wps_dir}/${metgrid_dir}/WRF/
                cp /home/${user}/climatex/remapping/compute_hourly.py ${wps_dir}/${metgrid_dir}/WRF/
                sed -i -e "s/JOBNAME/${user:0:2}H${ssyear:2}${ssmonth}_${eeyear:2}${eemonth}/g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_regrid_wrfhourly.sh
                sed -i -e "s%RUNDIR%${wps_dir}/${metgrid_dir}/WRF%g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_regrid_wrfhourly.sh
                sed -i -e "s/WPSDIR/${WPSDIR}/g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_regrid_wrfhourly.sh
                sed -i -e "s/METGRIDDIR/${metgrid_dir}/g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_regrid_wrfhourly.sh

                sbatch ${wps_dir}/${metgrid_dir}/WRF/submit_run_regrid_wrfhourly.sh
            fi
        fi
        
        ### REMAPPING OF WRFPRESS FILES ###
        squeue -u ${user} | grep ${user:0:2}P${ssyear:2}${ssmonth}_${eeyear:2}${eemonth}
        jobrunning=$?
        
        if [[ $jobrunning == 0 ]]; then 
            echo "Job ${user:0:2}P${ssyear:2}${ssmonth}_${eeyear:2}${eemonth} still running"
            echo ""  
        else
            if [ ! -e ${wps_dir}/${metgrid_dir}/WRF/RCLONE_WRFPRESS.OK ]; then
                cp /home/${user}/climatex/remapping/submit_run_regrid_wrfpress.sh ${wps_dir}/${metgrid_dir}/WRF/
                cp /home/${user}/climatex/remapping/run_regrid_wrfpress.sh ${wps_dir}/${metgrid_dir}/WRF/
                sed -i -e "s/JOBNAME/${user:0:2}P${ssyear:2}${ssmonth}_${eeyear:2}${eemonth}/g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_regrid_wrfpress.sh
                sed -i -e "s%RUNDIR%${wps_dir}/${metgrid_dir}/WRF%g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_regrid_wrfpress.sh
                sed -i -e "s/WPSDIR/${WPSDIR}/g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_regrid_wrfpress.sh
                sed -i -e "s/METGRIDDIR/${metgrid_dir}/g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_regrid_wrfpress.sh

                sbatch ${wps_dir}/${metgrid_dir}/WRF/submit_run_regrid_wrfpress.sh
            fi
        fi

    else

        # RCLONE FINISHED; skipping month, move on to the next one
        echo ""
        echo "${wps_dir}/${metgrid_dir}/WRF/RCLONE.OK found; moving on"
        echo ""
        prev_dir=${wps_dir}/${metgrid_dir}

        cd ${prev_dir}/WRF
      
        cd ${wps_dir}
        continue
    fi
done






