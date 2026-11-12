#!/bin/bash -l


module load cdo/2.0.5
SECONDS=0

user=`whoami`
slack=/home/${user}/climatex/slack.sh
username="$(whoami)@$(hostname -s)"
channel=climatex_runs
target='<!channel>'

if [ ! -e REMAPPED_nn ]; then mkdir REMAPPED_nn || echo hourly; fi
if [ ! -e REMAPPED_bic ]; then mkdir REMAPPED_bic || echo hourly; fi
if [ ! -e REMAPPED_bil ]; then mkdir REMAPPED_bil || echo hourly; fi
if [ ! -e REMAPPED_dis ]; then mkdir REMAPPED_dis || echo hourly; fi

REMAP_DIR=/home/${user}/climatex/remapping

doms=("d01" "d02" "d03")
#doms=("d01")

remap_files=`ls ${REMAP_DIR}/*remap.txt`

anycode=0

for each_dom in ${doms[@]}; do
    code=0
    mkdir -p REMAPPED_nn/${each_dom}/wrfhourly
    mkdir -p REMAPPED_bil/${each_dom}/wrfhourly
    mkdir -p REMAPPED_bic/${each_dom}/wrfhourly
    mkdir -p REMAPPED_dis/${each_dom}/wrfhourly

    date
    source ~/WRFV4.3.3/source_wrf
    time python3 ./compute_hourly.py ${each_dom} || code=1
    module load StdEnv/2020  intel/2020.1.217  openmpi/4.0.3 cdo/2.0.5
       
    for wrfhourly in `ls wrfhourly_${each_dom}*00`; do
 
        echo $wrfhourly

        for each_remap_txt in ${remap_files}; do
            stnname=`echo $each_remap_txt | rev |  cut -d'/' -f1 | rev |  cut -d'.' -f1`

            cdo -s remap,${each_remap_txt},${REMAP_DIR}/${stnname}_nn_${each_dom}.nc ${wrfhourly} REMAPPED_nn/${each_dom}/wrfhourly/${wrfhourly}.${stnname}.nn.nc || code=1 &
            sleep 0.25
            cdo -s remap,${each_remap_txt},${REMAP_DIR}/${stnname}_bil_${each_dom}.nc ${wrfhourly} REMAPPED_bil/${each_dom}/wrfhourly/${wrfhourly}.${stnname}.bil.nc || code=1 &
            sleep 0.25
            cdo -s remap,${each_remap_txt},${REMAP_DIR}/${stnname}_bic_${each_dom}.nc ${wrfhourly} REMAPPED_bic/${each_dom}/wrfhourly/${wrfhourly}.${stnname}.bic.nc || code=1 &
            sleep 0.25
            cdo -s remap,${each_remap_txt},${REMAP_DIR}/${stnname}_dis_${each_dom}.nc ${wrfhourly} REMAPPED_dis/${each_dom}/wrfhourly/${wrfhourly}.${stnname}.dis.nc || code=1 &
            sleep 0.25

            wait
        done
        
        if [[ $code == 1 ]]; then
            #$slack -c $channel -s HIGH -t "Problem with computing or interpolating hourly values" -u $username "`pwd`"
            anycode=1
        fi


    done
 
done


if [[ $anycode == 0 ]]; then
    touch WRFHOURLY.OK
else
    $slack -c $channel -s HIGH -t "Problem with computing or interpolating hourly values" -u $username "`pwd`"
fi

