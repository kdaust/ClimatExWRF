#!/bin/bash -l

for each in `ls -d /Volumes/LaCie/CLIMATEX/CMIP6/*/ | cut -d'/' -f6`; do 
    echo $each
    if [[ -e REGRIDDED/${each}.OK ]]; then
       echo "skipping $each"
    else
        if [[ $each == "GISS-E2-1-G" ]] || [[ $each == "INM-CM5-0" ]] || [[ $each == "IPSL-CM6A-LR" ]]; then
            echo "skipping $each cause something's wrong"
            listtos=`ls ./step_01_extract_deltas/Omon/${each}/*.nc`
            
            for eachtos in ${listtos[@]}; do
                echo $eachtos

                if [[ ! -e ${eachtos}_orig ]]; then
                    cp ${eachtos} ${eachtos}_orig
                fi
                
                cdo -selvar,tos ${eachtos}_orig temp.nc
                cdo remapnn,./step_01_extract_deltas/Omon/CanESM5/tos_CanESM5_pi_monthly_delta.nc temp.nc ${eachtos}
                rm -fv temp.nc
                
            done

            listsic=`ls ./step_01_extract_deltas/SImon/${each}/*.nc`
            
            for eachsic in ${listsic[@]}; do
                echo $eachsic

                if [[ ! -e ${eachsic}_orig ]]; then
                    cp ${eachsic} ${eachsic}_orig
                fi
                 
                cdo -selvar,siconc ${eachsic}_orig temp.nc
                cdo remapnn,./step_01_extract_deltas/SImon/CanESM5/siconc_CanESM5_pi_monthly_delta.nc temp.nc ${eachsic}
                rm -fv temp.nc
            done

            ./step_02_preproc_deltas.sh $each pi || exit 1
            touch REGRIDDED/${each}.OK

        else
            ./step_02_preproc_deltas.sh $each pi || exit 1
            touch REGRIDDED/${each}.OK
        fi
    fi
done
            

