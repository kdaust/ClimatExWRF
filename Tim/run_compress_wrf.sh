#!/bin/bash -l


module load cdo/2.0.5

user=`whoami`
slack=/home/${user}/climatex/slack.sh
username="$(whoami)@$(hostname -s)"
channel=climatex_runs
target='<!channel>'

rundir=`pwd`
anycode=0
SECONDS=0
for each_wrf in `ls wrfout*00`; do 
    code=0
    nccopy -d1 -s ./${each_wrf} ./${each_wrf}_compressed || code=1

    if [[ $code == 1 ]]; then
        echo "${each_wrf}_compressed failed to be produced"
        $slack -c $channel -s HIGH -t "compression failure in ${rundir}" -u $username "${rundir}/${each_wrf}_compressed unable to be produced during compression"
        anycode=1
    else
        echo "${each_wrf}_compressed produced"
        #rm -fv ./${each_wrf}
    fi
done

if [[ $anycode == 0 ]]; then
    touch COMPRESS.OK
fi



