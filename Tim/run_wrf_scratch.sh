#!/bin/bash -l
set -x
if [ $# -lt 6 ]; then
    # TODO: print usage
    echo "Usage: run_wrf SSYEAR SSMONTH SSDAY EEYEAR EEMONTH EEDAY"
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

get_restarts () {

    wrf_dir=${metgrid_dir}/WRF
    cd $wrf_dir
    
    d03_restart_n1=`ls wrfrst*d03* | tail -n1`
    d03_hour_n1=`echo $d03_restart_n1 | cut -d"_" -f4`
    d03_restart_n2=`ls wrfrst*d03* | tail -n2 | head -n1`

    d02_restart_n1=`ls wrfrst*d02* | tail -n1`
    d02_hour_n1=`echo $d02_restart_n1 | cut -d"_" -f4`
    d02_restart_n2=`ls wrfrst*d02* | tail -n2 | head -n1`

    d01_restart_n1=`ls wrfrst*d01* | tail -n1`
    d01_hour_n1=`echo $d01_restart_n1 | cut -d"_" -f4`
    d01_restart_n2=`ls wrfrst*d01* | tail -n2 | head -n1`

    if [ $d03_hour_n1 != $d02_hour_n1 ] || [ $d03_hour_n1 != $d01_hour_n1 ]; then
        d03_restart=$d03_restart_n2
        d02_restart=$d02_restart_n2
        d01_restart=$d01_restart_n2
    else
        d03_size_n1=`ls -l $d03_restart_n1 | cut -d' ' -f5`
        d03_size_n2=`ls -l $d03_restart_n2 | cut -d' ' -f5`

        d02_size_n1=`ls -l $d02_restart_n1 | cut -d' ' -f5`
        d02_size_n2=`ls -l $d02_restart_n2 | cut -d' ' -f5`

        d01_size_n1=`ls -l $d01_restart_n1 | cut -d' ' -f5`
        d01_size_n2=`ls -l $d01_restart_n2 | cut -d' ' -f5`

        if [ $d03_size_n1 != $d03_size_n2 ] || [ $d02_size_n1 != $d02_size_n2 ] || [ $d01_size_n1 != $d01_size_n2 ]; then
            d03_restart=$d03_restart_n2
            d02_restart=$d02_restart_n2
            d01_restart=$d01_restart_n2
        else
            d03_restart=$d03_restart_n1
            d02_restart=$d02_restart_n1
            d01_restart=$d01_restart_n1
        fi
    fi           
    restart_date=`echo $d03_restart | cut -d"_" -f3`
    restart_year=`echo $restart_date | cut -d"-" -f1`
    restart_month=`echo $restart_date | cut -d"-" -f2`
    restart_day=`echo $restart_date | cut -d"-" -f3`
    restart_hour=`echo $d03_restart | cut -d"_" -f4 | cut -d":" -f1`
    
    start_year=`grep "start_year" ${namelist_wrf}`
    start_month=`grep "start_month" ${namelist_wrf}`
    start_day=`grep "start_day" ${namelist_wrf}`
    start_hour=`grep "start_hour" ${namelist_wrf}`

    sed -i -e "s/${start_year}/start_year = ${restart_year},${restart_year},${restart_year},${restart_year},/g" ${namelist_wrf}
    sed -i -e "s/${start_month}/start_month = ${restart_month},${restart_month},${restart_month},${restart_month},/g" ${namelist_wrf}
    sed -i -e "s/${start_day}/start_day = ${restart_day},${restart_day},${restart_day},${restart_day},/g" ${namelist_wrf}
    sed -i -e "s/${start_hour}/start_hour = ${restart_hour},${restart_hour},${restart_hour},${restart_hour},/g" ${namelist_wrf}

    date_restart=$(date -u '+%s' -d "${restart_year}-${restart_month}-${restart_day} ${restart_hour}")
    date_end=$(date -u '+%s' -d "${eeyear}-${eemonth}-${eeday} ${eehour}")
    date_diff=`echo $(( ( date_end - date_restart )/(60*60) ))`

    sed -i -e "s/HHOURS/${date_diff}/g" ${namelist_wrf}

    if [[ 10#$restart_day > 02 ]]; then
        prev_restart_day=`printf %02d $(( 10#$restart_day - 1 ))`
        # Clean up restart files
        for each_wrfrst in `ls wrfrst*`; do
            echo $each_wrfrst
            old_rst_day=`echo $each_wrfrst | cut -d'_' -f3 | cut -d'-' -f3`
            if [[ 10#$old_rst_day != 10#$restart_day ]] && [[ 10#$old_rst_day != 10#$prev_restart_day ]]; then

                old_rst_hour=`echo $each_wrfrst | cut -d'_' -f4 | cut -d':' -f1`

                # Delete all previous (non-recent) dates NOT divisible by 3
                # Following if will evaluate to non-zero number if not divisible by 3
                if  (( 10#$old_rst_day % 3 )); then

                    echo $old_rst_day
                    ls $each_wrfrst >> delete_wrfrst.txt
                    rm -fv $each_wrfrst
                else
                    # Only delete hour 12 off of dates divisble by 3
                    if [[ $old_rst_hour != 12 ]]; then
                        ls $each_wrfrst >> delete_wrfrst.txt
                        rm -fv $each_wrfrst
                    else
                        echo "Skipping $each_wrfrst; day divisible by 3 but hour is 12"
                    fi

                fi

            else
                echo "Skipping ${each_wrfrst}; too recent"

            fi

        done
    fi

    cd $wps_dir

}

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
    namelist_wrf="namelist.wrf_${RUNNAME}"
else
    wps_dir=${ROOT_DIR}/WPS_${SSYEAR}${SSMONTH}_${EEYEAR}${EEMONTH}_${user}
    namelist_wrf="namelist.wrf"
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
    mkdir -p $metgrid_dir/WRF 
   
  
    # If wrf is finished, move on
    if [ ! -e ${wps_dir}/${metgrid_dir}/WRF/WRF.OK ]; then
        
        squeue -u ${user} | grep ${RUNNAME}${user:0:1}W${ssyear:2}${ssmonth}_${eeyear:2}${eemonth}
        jobrunning=$?
        
        # Don't do anything if a job is still running
        if [[ $jobrunning == 0 ]]; then
            echo ""
            echo "Job ${RUNNAME}${user:0:1}W${ssyear:2}${ssmonth}_${eeyear:2}${eemonth} still running; exiting"
            echo ""
            exit 0
        fi

        # Otherwise check that REAL.OK is produced
        SECONDS=0
        warning1=43200
        warning2=86400

        echo_flag=1
        while [ ! -e ${wps_dir}/${metgrid_dir}/WRF/REAL.OK ]; do
            echo "Waiting for ${wps_dir}/${metgrid_dir}/WRF/REAL.OK"
            duration=$SECONDS
          


            if [[ $echo_flag == 1 ]]; then 
                if [[ $SECONDS -ge $warning2 ]]; then
                    $slack -c $channel -s MED -t "${wps_dir}/${metgrid_dir}/WRF/REAL.OK not found for ${warning2} s" -u $username "$target ${wps_dir}/${metgrid_dir}"
                    echo "${warning2} s elapsed"
                    echo_flag=0
                elif
                   [[ $SECONDS -ge $warning1 ]]; then
                    $slack -c $channel -s MED -t "${wps_dir}/${metgrid_dir}/WRF/REAL.OK not found for ${warning1} s" -u $username "$target ${wps_dir}/${metgrid_dir}"
                    echo "${warning1} s elapsed"
                fi
            fi
            echo $SECONDS
            sleep 60
        done

        echo ""
        echo "${ssyear}_${ssmonth}"
        echo ""
 
        # Clean up namelist.input
        rm -f ${wps_dir}/${metgrid_dir}/WRF/namelist.input
        
        cp /home/${user}/climatex/${namelist_wrf} ${wps_dir}/${metgrid_dir}/WRF/${namelist_wrf}
        cp /home/${user}/climatex/submit_run_wrf.sh ${wps_dir}/${metgrid_dir}/WRF/submit_run_wrf.sh
        
        ### Edit namelist.wrf ### 
        sed -i -e "s/SSYEAR/${ssyear}/g" ${wps_dir}/${metgrid_dir}/WRF/${namelist_wrf}
        sed -i -e "s/SSMONTH/${ssmonth}/g" ${wps_dir}/${metgrid_dir}/WRF/${namelist_wrf}
        #sed -i -e "s/EEYEAR/${eeyear}/g" ${wps_dir}/${metgrid_dir}/WRF/${namelist_wrf}
        
        # Handle situation if end day is not 01 (i.e. last dir finished in same month)
        if [[ ${metgrid_dir} == metgrid_${EEYEAR}_${EEMONTH} ]]; then
            #eeday=${EEDAY}
            #eemonth=${ssmonth}

            eeday=01

            if [[ $EEMONTH == 12 ]]; then
                eemonth=01
                eeyear=$(( EEYEAR + 1 ))
            else
                eemonth_num=$(( EEMONTH + 1 ))
                eemonth=`printf "%02d" $eemonth_num`
                eeyear=$EEYEAR
            fi

        else
            eeday=01
        fi
        
        sed -i -e "s/EEYEAR/${eeyear}/g" ${wps_dir}/${metgrid_dir}/WRF/${namelist_wrf}
        sed -i -e "s/EEMONTH/${eemonth}/g" ${wps_dir}/${metgrid_dir}/WRF/${namelist_wrf}
        sed -i -e "s/EEDAY/${eeday}/g" ${wps_dir}/${metgrid_dir}/WRF/${namelist_wrf}

        # Only model end in December finishes 12 hours into January
        if [[ ${metgrid_dir} == metgrid_${MYEAR}_12 ]]; then
            eehour=12
        else
            eehour=00
        fi

        sed -i -e "s/EEHOUR/${eehour}/g" ${wps_dir}/${metgrid_dir}/WRF/${namelist_wrf}

        if [[ ${metgrid_dir} == metgrid_${SSYEAR}_${SSMONTH} ]] && [ ! -e ${wps_dir}/${metgrid_dir}/WRF/STATIC.OK ]; then
            sed -i -e "s/SSDAY/01/g" ${wps_dir}/${metgrid_dir}/WRF/${namelist_wrf}
            sed -i -e "s/SSDAY/01/g" ${wps_dir}/${metgrid_dir}/WRF/${namelist_wrf}
            sed -i -e "s/SSDAY/01/g" ${wps_dir}/${metgrid_dir}/WRF/${namelist_wrf}
            sed -i -e "s/SSDAY/01/g" ${wps_dir}/${metgrid_dir}/WRF/${namelist_wrf}
               
            sed -i -e "s/SSHOUR1/00/g" ${wps_dir}/${metgrid_dir}/WRF/${namelist_wrf}
            sed -i -e "s/SSHOUR2/03/g" ${wps_dir}/${metgrid_dir}/WRF/${namelist_wrf}
            sed -i -e "s/SSHOUR3/06/g" ${wps_dir}/${metgrid_dir}/WRF/${namelist_wrf}
            sed -i -e "s/SSHOUR4/09/g" ${wps_dir}/${metgrid_dir}/WRF/${namelist_wrf}

            sed -i -e "s/HHOURS/12/g" ${wps_dir}/${metgrid_dir}/WRF/${namelist_wrf}
            sed -i -e "s/ADAPTIVE/.false./g" ${wps_dir}/${metgrid_dir}/WRF/${namelist_wrf}
            sed -i -e "s/RESTART/.false./g" ${wps_dir}/${metgrid_dir}/WRF/${namelist_wrf}
        else
            
            last_wrfrst_d01=`ls ${prev_dir}/WRF/wrfrst*d01* | tail -n1`
            cp -v  $last_wrfrst_d01 ${wps_dir}/${metgrid_dir}/WRF/
            last_wrfrst_d02=`ls ${prev_dir}/WRF/wrfrst*d02* | tail -n1`
            cp -v  $last_wrfrst_d02 ${wps_dir}/${metgrid_dir}/WRF/
            last_wrfrst_d03=`ls ${prev_dir}/WRF/wrfrst*d03* | tail -n1`
            cp -v  $last_wrfrst_d03 ${wps_dir}/${metgrid_dir}/WRF/
     
            # Edit namelist.input with latest valid wrfrst files
            get_restarts
            
            sed -i -e "s/ADAPTIVE/.true./g" ${wps_dir}/${metgrid_dir}/WRF/${namelist_wrf}
            sed -i -e "s/RESTART/.true./g" ${wps_dir}/${metgrid_dir}/WRF/${namelist_wrf}


        fi


        sed -i -e "s/JOBNAME/${RUNNAME}${user:0:1}W${ssyear:2}${ssmonth}_${eeyear:2}${eemonth}/g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_wrf.sh
        sed -i -e "s%RUNDIR%${wps_dir}/${metgrid_dir}/WRF%g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_wrf.sh
        
        squeue -u ${user} | grep ${RUNNAME}${user:0:1}W${ssyear:2}${ssmonth}_${eeyear:2}${eemonth}
        jobrunning=$?
        
        # Don't do anything if a job is still running
        if [[ $jobrunning == 0 ]]; then
            echo ""
            echo "Job ${RUNNAME}${user:0:1}W${ssyear:2}${ssmonth}_${eeyear:2}${eemonth} still running; exiting"
            echo ""
            exit 0
        fi
     

        cp ${wps_dir}/${metgrid_dir}/WRF/${namelist_wrf} ${wps_dir}/${metgrid_dir}/WRF/namelist.input
        
        nonamelist=0
        cat ${wps_dir}/${metgrid_dir}/WRF/namelist.input || nonamelist=1

        if [[ $nonamelist == 1 ]]; then
            echo "Something wrong with namelist.input; exiting"
            exit 1
        fi

        echo "Submitting job ${RUNNAME}${user:0:1}W${ssyear:2}${ssmonth}_${eeyear:2}${eemonth} to run wrf.exe"
        
        cd ${wps_dir}/${metgrid_dir}/WRF
        sbatch submit_run_wrf.sh
        cd ${wps_dir}    
        
        echo "Waiting for ${RUNNAME}${user:0:1}W${ssyear:2}${ssmonth}_${eeyear:2}${eemonth} to finish; exiting in the meantime"
        echo ""
        exit 0

    else
        echo ""
        echo "${wps_dir}/${metgrid_dir}/WRF/WRF.OK found; moving on"
        echo ""
        prev_dir=${wps_dir}/${metgrid_dir}

        cd ${prev_dir}/WRF


        cd ${wps_dir}
        continue
    fi
done






