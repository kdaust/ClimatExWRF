#!/bin/bash -l
#
#

rm -fv *inventory*.csv

while read line; do

    server=`echo $line | cut -d',' -f1`
    directory=`echo $line | cut -d',' -f2 | xargs`
    if [[ $server == climatex ]] || [[ $server == climatex2 ]]; then
        echo $directory 
        wy=`echo "${directory:11:4}"`
        #echo $wy
        wym1=`echo "${directory:4:4}"`
        #echo $wym1
        months=(${wym1}_09 ${wym1}_10 ${wym1}_11 ${wym1}_12 ${wy}_01 ${wy}_02 ${wy}_03 ${wy}_04 ${wy}_05 ${wy}_06 ${wy}_07 ${wy}_08 ${wy}_09)
        #echo ${months[@]}
        for month in ${months[@]}; do
            check=`rclone lsf ${server}:Share_Data/${directory}/metgrid_${month}/WRFOUT --include=wrf{out,rst}_d03*-01_00:00:00*`
            code=0
            if [[ -z $check ]]; then
                echo ${server}:Share_Data/${directory}/metgrid_${month}/WRFOUT missing
                #exit 1
                code=1
            fi
            echo ${server}:Share_Data/${directory}/metgrid_${month}/WRFOUT >> wrfout_inventory.csv
            echo ${server}:Share_Data/${directory}/metgrid_${month}/WRFOUT,$code >> wrfout_inventory_status.csv
        done
    else
        first=`echo $server | cut -d';' -f1`
        second=`echo $server | cut -d';' -f2`

        first_server=`echo $first | cut -d':' -f1 | xargs`
        first_months=`echo $first | cut -d':' -f2`

        second_server=`echo $second | cut -d':' -f1 | xargs`
        second_months=`echo $second | cut -d':' -f2`

        #echo $first_server $first_months
        #echo $second_server $second_months
        echo $directory
        for month in ${first_months[@]}; do
            check=`rclone lsf ${first_server}:Share_Data/${directory}/metgrid_${month}/WRFOUT --include=wrfout_d03*-01_00:00:00*compressed`
            code=0
            if [[ -z $check ]]; then
                echo ${first_server}:Share_Data/${directory}/metgrid_${month}/WRFOUT missing
                #exit 1
                code=1
            fi
            echo ${first_server}:Share_Data/${directory}/metgrid_${month}/WRFOUT >> wrfout_inventory.csv
            echo ${first_server}:Share_Data/${directory}/metgrid_${month}/WRFOUT,$code >> wrfout_inventory_status.csv
        done

        for month in ${second_months[@]}; do
            check=`rclone lsf ${second_server}:Share_Data/${directory}/metgrid_${month}/WRFOUT --include=wrfout_d03*-01_00:00:00*compressed`
            code=0
            if [[ -z $check ]]; then
                echo ${second_server}:Share_Data/${directory}/metgrid_${month}/WRFOUT missing
                #exit 1
                code=1
            fi
            echo ${second_server}:Share_Data/${directory}/metgrid_${month}/WRFOUT >> wrfout_inventory.csv
            echo ${second_server}:Share_Data/${directory}/metgrid_${month}/WRFOUT,$code >> wrfout_inventory_status.csv
        done

    fi

done < test_directories.csv
