#!/bin/bash
#SBATCH --time=22:59:00
#SBATCH --ntasks-per-node=192
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=3950M
#SBATCH --job-name=JOBNAME
#SBATCH --account=rrg-rstull

ulimit -s unlimited
ulimit -v unlimited
ulimit -m unlimited
module load StdEnv/2023  gcc/13.3 openmpi/5.0.3

get_restarts () {

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
    
    start_year=`grep "start_year" namelist.input`
    start_month=`grep "start_month" namelist.input`
    start_day=`grep "start_day" namelist.input`
    start_hour=`grep "start_hour" namelist.input`
    hhours=`grep "run_hours" namelist.input`

    sed -i -e "s/${start_year}/start_year = ${restart_year},${restart_year},${restart_year},${restart_year},/g" namelist.input
    sed -i -e "s/${start_month}/start_month = ${restart_month},${restart_month},${restart_month},${restart_month},/g" namelist.input
    sed -i -e "s/${start_day}/start_day = ${restart_day},${restart_day},${restart_day},${restart_day},/g" namelist.input
    sed -i -e "s/${start_hour}/start_hour = ${restart_hour},${restart_hour},${restart_hour},${restart_hour},/g" namelist.input

    date_restart=$(date -u '+%s' -d "${restart_year}-${restart_month}-${restart_day} ${restart_hour}")
    
    eeyear=`grep "end_year" namelist.input | cut -d'=' -f2 | cut -d',' -f1 | xargs`
    eemonth=`grep "end_month" namelist.input | cut -d'=' -f2 | cut -d',' -f1 | xargs`
    eeday=`grep "end_day" namelist.input | cut -d'=' -f2 | cut -d',' -f1 | xargs`
    eehour=`grep "end_hour" namelist.input | cut -d'=' -f2 | cut -d',' -f1 | xargs`

    date_end=$(date -u '+%s' -d "${eeyear}-${eemonth}-${eeday} ${eehour}")
    date_diff=`echo $(( ( date_end - date_restart )/(60*60) ))`

    sed -i -e "s/${hhours}/run_hours = ${date_diff}/g" namelist.input

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
                    #rm -fv $each_wrfrst

                else
                    # Only delete hour 12 off of dates divisble by 3
                    if [[ $old_rst_hour != 12 ]]; then
                        ls $each_wrfrst >> delete_wrfrst.txt
                        #rm -fv $each_wrfrst
                    else
                        echo "Skipping $each_wrfrst; day divisible by 3 but hour is 12"
                    fi

                fi

            else
                echo "Skipping ${each_wrfrst}; too recent"
            
            fi

        done
    fi

}

user=`whoami`
slack=/home/${user}/climatex/slack.sh
username="$(whoami)@$(hostname -s)"
channel=climatex_runs
target='<!channel>'



#source /home/${user}/WRFV4.3.3/source_wrf
rundir=RUNDIR

cd $rundir

last_wrf=`ls wrfout*d03* | tail -n1`

$slack -c $channel -s MED -t "JOBNAME, LAST_WRF=${last_wrf}" -u $username "$target RUNDIR"


code=0

SECONDS=0
export I_MPI_TUNING_BIN=${rundir}/tuning.dat
time srun ${rundir}/wrf.exe || code=1
duration=$SECONDS

touch ${duration}.OK

if [ $code == 1 ]; then
    
    count=0
    
    while [ $count -le 6 ]; do
        wrfout=`ls wrfout* | tail -n1`
        $slack -c $channel -s MED -t "JOBNAME has failed to complete wrf.exe after ${duration} s; reducing CFL" -u $username "$target Last wrfout: ${rundir}/$wrfout"
        
        target_cfl=`grep target_cfl namelist.input`
        if [ $count == 0 ]; then
            new_target_cfl="target_cfl=1.05,1.15,1.25,1.3,"
        elif [ $count == 1 ]; then
            new_target_cfl="target_cfl=1.0,1.1,1.2,1.25,"
        elif [ $count == 2 ]; then
            new_target_cfl="target_cfl=0.95,1.05,1.15,1.2,"
        elif [ $count == 3 ]; then
            new_target_cfl="target_cfl=0.9,1.0,1.1,1.15,"  
        elif [ $count == 4 ]; then
            new_target_cfl="target_cfl=0.85,0.95,1.05,1.1,"
        elif [ $count == 5 ]; then
            new_target_cfl="target_cfl=0.8,0.9,1.0,1.05,"
        else
            new_target_cfl="target_cfl=0.7,0.8,0.9,0.95,"
        fi
      
        sed -i -e  "s/${target_cfl}/${new_target_cfl}/g" namelist.input
        
        restart_line=`grep "restart" namelist.input | head -n1`
        sed -i -e "s/${restart_line}/restart = .true.,/g" namelist.input
        
        # Edit namelist.input with latest valid wrfrst files
        get_restarts
        
        code=0
     
        SECONDS=0
        time srun ${rundir}/wrf.exe || code=1
        duration=$SECONDS

        if [ $code == 0 ]; then
            $slack -c $channel -s LOW -t "JOBNAME COMPLETED in ${duration} s after $count additional retries" -u $username "$target $rundir"
            touch SUCCESS_${duration}.OK
            touch WRF.OK
            exit 0  
        fi
        touch RESTARTED_0${count}_${duration}.OK
        count=$(( count + 1 ))

    done
    wrfout=`ls wrfout* | tail -n1`
    $slack -c $channel -s HIGH -t "JOBNAME has FAILED to complete wrf.exe after ${duration} s" -u $username "$target Last wrfout: ${rundir}/$wrfout"
    touch FAIL_${duration}.OK
    exit 1
else

    $slack -c $channel -s LOW -t "JOBNAME COMPLETED in ${duration} s" -u $username "$target $rundir"
    touch SUCCESS_${duration}.OK
    is_adaptive=`grep "use_adaptive_time_step" namelist.input | cut -d'=' -f2`
    if [[  $is_adaptive == *".false.,"* ]]; then
        touch STATIC.OK
    else
        touch WRF.OK
    fi
    exit 0
fi 


