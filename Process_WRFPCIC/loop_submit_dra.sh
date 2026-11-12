#!/bin/bash -l
#

lom=(metgrid_2020_11 metgrid_2020_12 metgrid_2021_01 metgrid_2021_02 metgrid_2021_03 metgrid_2021_04 metgrid_2021_05 metgrid_2021_06 metgrid_2021_07 metgrid_2021_08 metgrid_2021_09)

root_dir=/home/tchui/scratch/PCIC/MAKE_PCIC/
scripts_dir=/home/tchui/scratch/PCIC/ClimatExWRF/Process_WRFPCIC

for each_metgrid in ${lom[@]}; do
    echo ""
    echo "=================="
    echo $each_metgrid
    echo ""
    mkdir -p ${root_dir}/${each_metgrid}
    ls -lh ${root_dir}/${each_metgrid}
    cd ${root_dir}/${each_metgrid}
    cp -v ${scripts_dir}/*py .
    cp -v ${scripts_dir}/*csv .

    cp -v ${scripts_dir}/process_pcic.sh .
    cp -v ${scripts_dir}/submit_dra.sh .
    
    sed -i s/JOBNAME/${each_metgrid}/g submit_dra.sh
    chmod 777 *
    sbatch submit_dra.sh
    cd ${root_dir}
    echo ""
done
