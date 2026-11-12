#!/bin/bash -l

if [ $# -lt 6 ]; then
    # TODO: print usage
    echo "Usage: run_real SSYEAR SSMONTH SSDAY EEYEAR EEMONTH EEDAY"
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
    namelist_real="namelist.real_${RUNNAME}"
    WPSDIR=WPS_${SSYEAR}${SSMONTH}_${EEYEAR}${EEMONTH}_${user}_${RUNNAME}
else
    wps_dir=${ROOT_DIR}/WPS_${SSYEAR}${SSMONTH}_${EEYEAR}${EEMONTH}_${user}
    namelist_real="namelist.real"
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
    mkdir -p $metgrid_dir/WRF 
   
  
    # If real is finished, skip this process
    if [ ! -e ${metgrid_dir}/WRF/REAL.OK ]; then
      
        SECONDS=0
        warning1=43200
        warning2=86400

        echo_flag=1
        while [ ! -e ${metgrid_dir}/METGRID.OK ]; do
            echo "Waiting for ${metgrid_dir}/METGRID.OK"
            duration=$SECONDS
            
            if [[ $SECONDS -ge $warning2 ]] && [[ $echo_flag == 1 ]]; then
                $slack -c $channel -s MED -t "${metgrid_dir}/METGRID.OK not found for ${warning2} s" -u $username "$target ${wps_dir}/${metgrid_dir}"
                echo "${warning2} s elapsed"
                echo_flag=0
            elif
               [[ $SECONDS -ge $warning1 ]] && [[ $echo_flag == 1 ]]; then
                $slack -c $channel -s MED -t "${metgrid_dir}/METGRID.OK not found for ${warning1} s" -u $username "$target ${wps_dir}/${metgrid_dir}"
                echo "${warning1} s elapsed"
            fi
            echo $SECONDS
            sleep 60
        done

        echo ""
        echo "${ssyear}_${ssmonth}"
        echo "Populating ${wps_dir}/${metgrid_dir}/WRF"
        echo ""

        ln -sf /home/${user}/WRFV4.3.3/WRF-4.4.2/test/em_real/* ${metgrid_dir}/WRF/
        if [[ $RUNNAME == pi ]]; then
            rm -fv ${metgrid_dir}/WRF/CAMtr_volume_mixing_ratio
            cp -v /home/${user}/WRFV4.3.3/WRF-4.4.2/test/em_real/CAMtr_volume_mixing_ratio.SSP245_pi ${metgrid_dir}/WRF/CAMtr_volume_mixing_ratio
        elif [[ $RUNNAME == 2 ]]; then
            rm -fv ${metgrid_dir}/WRF/CAMtr_volume_mixing_ratio
            cp -v /home/${user}/WRFV4.3.3/WRF-4.4.2/test/em_real/CAMtr_volume_mixing_ratio.SSP585_2 ${metgrid_dir}/WRF/CAMtr_volume_mixing_ratio
        elif [[ $RUNNAME == 3 ]]; then
            rm -fv ${metgrid_dir}/WRF/CAMtr_volume_mixing_ratio
            cp -v /home/${user}/WRFV4.3.3/WRF-4.4.2/test/em_real/CAMtr_volume_mixing_ratio.SSP585_3 ${metgrid_dir}/WRF/CAMtr_volume_mixing_ratio
        fi


        rm -f ${metgrid_dir}/WRF/${namelist_real}
        rm -f ${metgrid_dir}/WRF/namelist.input

        cp /home/${user}/climatex/${namelist_real} ${metgrid_dir}/WRF/${namelist_real}
        cp /home/${user}/climatex/submit_run_real.sh ${metgrid_dir}/WRF/submit_run_real.sh

        sed -i -e "s/SSYEAR/${ssyear}/g" ${metgrid_dir}/WRF/${namelist_real}
        sed -i -e "s/SSMONTH/${ssmonth}/g" ${metgrid_dir}/WRF/${namelist_real}
        sed -i -e "s/SSDAY/01/g" ${metgrid_dir}/WRF/${namelist_real}
        sed -i -e "s/SSDAY/01/g" ${metgrid_dir}/WRF/${namelist_real}
        sed -i -e "s/SSDAY/01/g" ${metgrid_dir}/WRF/${namelist_real}
        sed -i -e "s/SSDAY/01/g" ${metgrid_dir}/WRF/${namelist_real}
        
        # Only model start has a staggered start for nests and static time step
        if [[ ${metgrid_dir} == metgrid_${SSYEAR}_${SSMONTH} ]]; then
            sed -i -e "s/SSHOUR1/00/g" ${metgrid_dir}/WRF/${namelist_real}
            sed -i -e "s/SSHOUR2/03/g" ${metgrid_dir}/WRF/${namelist_real}
            sed -i -e "s/SSHOUR3/06/g" ${metgrid_dir}/WRF/${namelist_real}
            sed -i -e "s/SSHOUR4/09/g" ${metgrid_dir}/WRF/${namelist_real}
            
            sed -i -e "s/ADAPTIVE/.false./g" ${metgrid_dir}/WRF/${namelist_real}
            sed -i -e "s/RESTART/.false./g" ${metgrid_dir}/WRF/${namelist_real}

        else
            sed -i -e "s/SSHOUR1/00/g" ${metgrid_dir}/WRF/${namelist_real}
            sed -i -e "s/SSHOUR2/00/g" ${metgrid_dir}/WRF/${namelist_real}
            sed -i -e "s/SSHOUR3/00/g" ${metgrid_dir}/WRF/${namelist_real}
            sed -i -e "s/SSHOUR4/00/g" ${metgrid_dir}/WRF/${namelist_real}

            sed -i -e "s/ADAPTIVE/.true./g" ${metgrid_dir}/WRF/${namelist_real}
            sed -i -e "s/RESTART/.true./g" ${metgrid_dir}/WRF/${namelist_real}

        fi

        sed -i -e "s/EEYEAR/${eeyear}/g" ${metgrid_dir}/WRF/${namelist_real}
        
        # Handle situation if end day is not 01 (i.e. last dir finished in same month)
        if [[ ${metgrid_dir} == metgrid_${EEYEAR}_${EEMONTH} ]]; then
            eeday=${EEDAY}
            eemonth=${ssmonth}
        else
            eeday=01
        fi
        
        sed -i -e "s/EEMONTH/${eemonth}/g" ${metgrid_dir}/WRF/${namelist_real}
        sed -i -e "s/EEDAY/${eeday}/g" ${metgrid_dir}/WRF/${namelist_real}

        # Only model end in December finishes 12 hours into January
        if [[ ${metgrid_dir} == metgrid_${MYEAR}_12 ]]; then
            eehour=12
        else
            eehour=00
        fi

        sed -i -e "s/EEHOUR/${eehour}/g" ${metgrid_dir}/WRF/${namelist_real}

        date_start=$(date -u '+%s' -d "${ssyear}-${ssmonth}-01 00")
        date_end=$(date -u '+%s' -d "${eeyear}-${eemonth}-${eeday} ${eehour}")
        date_diff=`echo $(( ( date_end - date_start )/(60*60) ))`
        
        echo "${ssyear}-${ssmonth}-01 00"
        echo "${eeyear}-${eemonth}-${eeday} ${eehour}"
        echo $date_start
        echo $date_end
        echo $date_diff
        
        sed -i -e "s/HHOURS/${date_diff}/g" ${metgrid_dir}/WRF/${namelist_real}
       

        sed -i -e "s/JOBNAME/${user:0:2}R${ssyear:2}${ssmonth}_${eeyear:2}${eemonth}/g" ${metgrid_dir}/WRF/submit_run_real.sh
        sed -i -e "s%RUNDIR%${wps_dir}/${metgrid_dir}/WRF%g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_real.sh
        sed -i -e "s/WPSDIR/${WPSDIR}/g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_real.sh
        sed -i -e "s/METGRIDDIR/${metgrid_dir}/g" ${wps_dir}/${metgrid_dir}/WRF/submit_run_real.sh
 
        # Clean up old met_em links, and re-link
        rm -f ${metgrid_dir}/WRF/met_em*d0*
        ln -s ${wps_dir}/${metgrid_dir}/met_em*d0* ${metgrid_dir}/WRF/
      

        cp ${metgrid_dir}/WRF/${namelist_real} ${metgrid_dir}/WRF/namelist.input
        echo "Submitting job ${user:0:2}R${ssyear:2}${ssmonth}_${eeyear:2}${eemonth} to run real.exe"
        
        cd ${metgrid_dir}/WRF
        sbatch submit_run_real.sh
        cd ${wps_dir}    
          
    else
        echo ""
        echo "${metgrid_dir}/WRF/REAL.OK found; skipping directory creation and running real"
        echo ""
    fi
done






