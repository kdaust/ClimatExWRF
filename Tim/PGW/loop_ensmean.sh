#!/bin/bash -l

variables=( tos tas hurs ps ta hur ua va zg siconc ts )
#variables=( ta )
for each_var in ${variables[@]}; do

    varlist=()
    for each in `ls -d /Volumes/LaCie/CLIMATEX/CMIP6/*/ | cut -d'/' -f6`; do 

        varfile=`ls REGRIDDED/${each}_pi/${each_var}_*delta.nc`
        varlist+=("$varfile")

    done
    echo ${varlist[@]}
    cdo ensmean ${varlist[@]} REGRIDDED/ENS_pi/${each_var}_ENS_pi_monthly_delta.nc 
    
    varlist=()
    for each in `ls -d /Volumes/LaCie/CLIMATEX/CMIP6/*/ | cut -d'/' -f6`; do 

        varfile=`ls REGRIDDED/${each}_pi/${each_var}_*deg1*.nc`
        varlist+=("$varfile")

    done
    echo ${varlist[@]}
    cdo ensmean ${varlist[@]} REGRIDDED/ENS_pi/${each_var}_ENS_pi_deg1_monthly.nc 

done                

