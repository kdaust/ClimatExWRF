#!/bin/bash -l


module load cdo/2.0.5

user=`whoami`
slack=/home/${user}/climatex/slack.sh
username="$(whoami)@$(hostname -s)"
channel=climatex_runs
target='<!channel>'


SECONDS=0

if [ ! -e REMAPPED_nn ]; then mkdir REMAPPED_nn || echo hello; fi
if [ ! -e REMAPPED_bic ]; then mkdir REMAPPED_bic || echo hello; fi
if [ ! -e REMAPPED_bil ]; then mkdir REMAPPED_bil || echo hello; fi
if [ ! -e REMAPPED_dis ]; then mkdir REMAPPED_dis || echo hello; fi

REMAP_DIR=/home/${user}/climatex/remapping

doms=("d01" "d02" "d03")
#doms=("d01")

remap_files=`ls ${REMAP_DIR}/*remap.txt`
anycode=0
for each_dom in ${doms[@]}; do
    mkdir -p REMAPPED_nn/${each_dom}/wrfuvic
    mkdir -p REMAPPED_bil/${each_dom}/wrfuvic
    mkdir -p REMAPPED_bic/${each_dom}/wrfuvic
    mkdir -p REMAPPED_dis/${each_dom}/wrfuvic

    date
    
    for wrfuvic in `ls wrfuvic_${each_dom}*00`; do
 
        echo $wrfuvic

        for each_remap_txt in ${remap_files}; do
            stnname=`echo $each_remap_txt | rev |  cut -d'/' -f1 | rev |  cut -d'.' -f1`

            cdo -s remap,${each_remap_txt},${REMAP_DIR}/${stnname}_nn_${each_dom}.nc ${wrfuvic} REMAPPED_nn/${each_dom}/wrfuvic/${wrfuvic}.${stnname}.nn.nc || code=1 &
            sleep 0.25
            cdo -s remap,${each_remap_txt},${REMAP_DIR}/${stnname}_bil_${each_dom}.nc ${wrfuvic} REMAPPED_bil/${each_dom}/wrfuvic/${wrfuvic}.${stnname}.bil.nc || code=1 &
            sleep 0.25
            cdo -s remap,${each_remap_txt},${REMAP_DIR}/${stnname}_bic_${each_dom}.nc ${wrfuvic} REMAPPED_bic/${each_dom}/wrfuvic/${wrfuvic}.${stnname}.bic.nc || code=1 &
            sleep 0.25
            cdo -s remap,${each_remap_txt},${REMAP_DIR}/${stnname}_dis_${each_dom}.nc ${wrfuvic} REMAPPED_dis/${each_dom}/wrfuvic/${wrfuvic}.${stnname}.dis.nc || code=1 &
            sleep 0.25


            wait
        done

        if [[ $code == 1 ]]; then
            $slack -c $channel -s HIGH -t "Problem with interpolating uvic values" -u $username "`pwd`"
            anycode=1
        fi

    done

    #comp=0
    #tar -czvf REMAPPED_nn/${each_dom}.tar.gz REMAPPED_nn/${each_dom} || comp=1

    #if [[ $comp == 1 ]]; then
    #    echo "ERROR WITH COMPRESSING REMAPPED_nn/${each_dom}"
    #else
    #    rm -rfv REMAPPED_nn/${each_dom}
    #fi

done


if [[ $anycode == 0 ]]; then
    touch WRFUVIC.OK
fi
