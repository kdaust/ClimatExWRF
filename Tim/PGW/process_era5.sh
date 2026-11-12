#!/bin/bash -l

module load cdo/2.0.5

#directories=(2020_09 2020_10 2020_11 2020_12 2021_01 2021_02 2021_03 2021_04 2021_05 2021_07 2021_08 2021_09)
#directory=(2021_04)
indir=$1
directories=("$indir")

for each_dir in ${directories[@]}; do 
    cd $each_dir

    for each in `ls ERA5*single*grib`; do 
        yyyy=`echo $each | cut -d'_' -f3`
        mm=`echo $each | cut -d'_' -f4`
        dd=`echo $each | cut -d'_' -f5`
        hh=`echo $each | cut -d'_' -f6 | cut -d':' -f1`
        each_pres=`ls ERA5*pressure*${yyyy}_${mm}_${dd}_${hh}*grib`
        echo $each_pres
        datestr=${yyyy}${mm}${dd}${hh}
        
        rm -fv temp.grib
        cdo merge $each $each_pres temp.grib
        cdo -f nc copy temp.grib ERA5_final_${datestr}.nc

    done

    cd ..
    ./step_03_apply_to_era_pres.sh ENS pi ${yyyy} ${mm}

    cd FINAL/ENS_pi

    module load cdo/2.0.5
    for each in `ls ERA5_final_${yyyy}${mm}*nc`; do cdo -f grb1 copy $each ${each}.grib; done
done



