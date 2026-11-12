#!/bin/bash -l

module load StdEnv/2020  intel/2020.1.217  openmpi/4.0.3
module load cdo/2.0.5

user=`whoami`
slack=/home/${user}/climatex/slack.sh
username="$(whoami)@$(hostname -s)"
channel=climatex_runs
target='<!channel>'


SECONDS=0


REMAP_DIR=/home/${user}/climatex/remapping

doms=("d01" "d02" "d03")
anycode=0
for each_dom in ${doms[@]}; do

    date
    
    for wrfpress in `ls wrfpress_${each_dom}*00`; do
        
        code=0
        echo $wrfpress
        cdo -sellevel,1 -selname,GHT_PL,U_PL,V_PL,T_PL,TD_PL,Q_PL,RH_PL ./${wrfpress} ./GHT85_${wrfpress}.nc || code=1
        cdo -sellevel,3 -selname,GHT_PL,U_PL,V_PL,T_PL,TD_PL,Q_PL,RH_PL ./${wrfpress} ./GHT50_${wrfpress}.nc || code=1
        cdo -sellevel,5 -selname,GHT_PL,U_PL,V_PL,T_PL,TD_PL,Q_PL,RH_PL ./${wrfpress} ./GHT25_${wrfpress}.nc || code=1


        #if [[ $code == 1 ]]; then
        #    #$slack -c $channel -s HIGH -t "Problem with interpolating wrfpress values" -u $username "`pwd`"
        #    anycode=1
        #fi

      
    done


done


if [[ $anycode == 0 ]]; then
    touch WRFPRESS.OK
fi
