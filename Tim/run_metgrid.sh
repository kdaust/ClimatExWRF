#!/bin/bash -l
set -x 
if [ $# -lt 6 ]; then
    # TODO: print usage
    echo "Usage: run_metgrid SSYEAR SSMONTH SSDAY EEYEAR EEMONTH EEDAY"
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

endyearp1=$(( EEYEAR + 1 ))

ROOT_DIR=/scratch/${user}/CLIMATEX

if [ ! -z "$RUNNAME" ]; then
    ungrib_dir=${ROOT_DIR}/WPS_${SSYEAR}${SSMONTH}_${EEYEAR}${EEMONTH}_${user}_${RUNNAME}
else
    ungrib_dir=${ROOT_DIR}/WPS_${SSYEAR}${SSMONTH}_${EEYEAR}${EEMONTH}_${user}
fi

IBCS_DIR=${ungrib_dir}/ibcs/

cd ${ungrib_dir} || exit 1

SECONDS=0
warning1=21600
warning2=43200
echo_flag=1
while [ ! -e UNGRIB.OK ]; do
    echo "Waiting for ${ungrib_dir}/UNGRIB.OK"
    duration=$SECONDS
    
    if [[ $SECONDS -ge $warning2 ]] && [[ $echo_flag == 1 ]]; then
        $slack -c $channel -s MED -t "${ungrib_dir}/UNGRIB.OK not found for ${warning2} s" -u $username "$target ${ungrib_dir}"
        echo "${warning2} s elapsed"
        echo_flag=0
    elif
       [[ $SECONDS -ge $warning1 ]] && [[ $echo_flag == 1 ]]; then
        $slack -c $channel -s MED -t "${ungrib_dir}/UNGRIB.OK not found for ${warning1} s" -u $username "$target ${ungrib_dir}"
        echo "${warning1} s elapsed"
    fi
    echo $SECONDS
    sleep 60
done

curr_month=" "
for each_FILE in `ls FILE*`; do
    
    month=`echo $each_FILE | cut -d'-' -f2`
    year=`echo $each_FILE | cut -d'-' -f1 | cut -d':' -f2`
    day=`echo $each_FILE | cut -d'-' -f3 | cut -d'_' -f1`
    if [[ $month != $curr_month ]]; then
        curr_month=$month
        next_month=`printf %02d $(( 10#${curr_month} + 1 ))`

        if [[ $next_month == 13 ]]; then
            next_month=01
            endyear=$(( year + 1 ))

        else
            endyear=${year}
        fi

        mkdir metgrid_${year}_${curr_month}
        cd metgrid_${year}_${curr_month}

        if [ ! -e METGRID.OK ]; then
            
            rm -fv *met_em*

            ln -s /home/${user}/WRFV4.3.3/WPS-4.3.1/metgrid/metgrid.exe .
            ln -s ../{FILE,LAND,SNOW}*${year}-${curr_month}* .
            #ln -s ../{FILE,LAND,SNOW}*${EEYEAR}-01* .

            ln -s ../TAVGSFC .
            
            cp /home/${user}/WRFV4.3.3/WPS-4.3.1/metgrid/METGRID.TBL .
            cp /home/${user}/climatex/geo_em*.nc .
            if [[ $month == $SSMONTH ]] && [[ $year == $SSYEAR ]] && [[ -z $RUNNAME ]]; then 
                submit_metgrid_script=submit_run_metgrid_land.sh
            else
                submit_metgrid_script=submit_run_metgrid.sh
            fi

            cp /home/${user}/climatex/${submit_metgrid_script} .

            cp /home/${user}/climatex/namelist.metgrid ./namelist.metgrid
    
            if [[ $month == $SSMONTH ]] && [[ $year == $SSYEAR ]] && [[ -z $RUNNAME ]]; then
                cp /home/${user}/climatex/namelist.metgrid_land ./namelist.metgrid_land
                cp namelist.metgrid_land ./namelist.wps
                sed -i -e "s/SSYEAR/${year}/g" namelist.wps
                sed -i -e "s/SSMONTH/${curr_month}/g" namelist.wps
                sed -i -e "s/SSDAY/01/g" namelist.wps
                sed -i -e "s/EEYEAR/${year}/g" namelist.wps
                sed -i -e "s/EEMONTH/${curr_month}/g" namelist.wps
                sed -i -e "s/EEDAY/02/g" namelist.wps

                ./metgrid.exe || exit 1

                for each_met in `ls met_em*`; do
                    mv -v ./${each_met} ./land_${each_met}
                done

            fi


            cp namelist.metgrid namelist.wps
            sed -i -e "s/SSYEAR/${year}/g" namelist.wps
            sed -i -e "s/SSMONTH/${curr_month}/g" namelist.wps
            sed -i -e "s/SSDAY/01/g" namelist.wps
            sed -i -e "s/EEYEAR/${endyear}/g" namelist.wps
            

            if [[ $month == 12 ]] && [[ $endyear == $EEYEAR ]]; then
                sed -i -e "s/EEMONTH/${next_month}/g" namelist.wps
                sed -i -e "s/EEDAY/01/g" namelist.wps 
                ln -s ../FILE*${endyear}-01* .

            else
                sed -i -e "s/EEMONTH/${next_month}/g" namelist.wps
                sed -i -e "s/EEDAY/01/g" namelist.wps
                ln -s ../FILE*${year}*${next_month}*-01* .

            fi

            sed -i -e "s/METJOB/M${year}${curr_month}_${endyear}${next_month}/g" ${submit_metgrid_script}
            METDIR=`pwd`
            sed -i -e "s%METDIR%${METDIR}%g" ${submit_metgrid_script}

            if [[ $month == $SSMONTH ]] && [[ $year == $SSYEAR ]]; then
                sed -i -e "s%SSYEAR%${year}%g" ${submit_metgrid_script}
                sed -i -e "s%SSMONTH%${month}%g" ${submit_metgrid_script}
            fi 

            sbatch ${submit_metgrid_script}

        fi

        cd ..

    fi

done






