#!/bin/bash
#SBATCH --time=21:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=3950M
#SBATCH --job-name=JOBNAME
#SBATCH --account=def-rstull

ulimit -s unlimited
ulimit -v unlimited
ulimit -m unlimited

user=`whoami`
slack=/home/${user}/climatex/slack.sh
username="$(whoami)@$(hostname -s)"
channel=climatex_runs
target='<!channel>'

SECONDS=0
rundir=RUNDIR

$slack -c $channel -s MED -t "${rundir} started wrfout compression" -u $username "JOBNAME"

cd $rundir

if [[ ! -e COMPRESS.OK ]]; then
    time ./run_compress_wrf.sh >& COMPRESS.log
fi

if [[ ! -e COMPRESS.OK ]]; then
    exit 1 
fi


rclone mkdir climatex:Share_Data/WPSDIR/METGRIDDIR/WRFOUT
rclone copy -vc ./namelist.input climatex:Share_Data/WPSDIR/METGRIDDIR/WRFOUT/
for each_wrfinput in `wrf*_d0{1,2,3}`; do
    # Copy wrfinputs, wrfbdy and wrfdda to nextcloud
    rclone copy -vc ./${each_wrfinput} climatex:Share_Data/WPSDIR/METGRIDDIR/WRFOUT/
    sleep 1
done

for each_wrf in `ls ./wrfout*00`; do
    code=0
    echo $each_wrf

    # Copy to nextcloud
    rclone copy -vc ./${each_wrf}_compressed climatex:Share_Data/WPSDIR/METGRIDDIR/WRFOUT/
    sleep 1
    # Repeat to confirm no change
    rclone copy -vc ./${each_wrf}_compressed climatex:Share_Data/WPSDIR/METGRIDDIR/WRFOUT/

    # Check that it's properly on the remote server, and delete if it's there
    rclone ls climatex:Share_Data/WPSDIR/METGRIDDIR/WRFOUT/${each_wrf}_compressed || code=1

    if [[ $code == 0 ]]; then
        rm -v $each_wrf
    fi
done

for each_wrf in `ls ./wrfrst*_{00,12}*00`; do
    code=0
    echo $each_wrf

    #nccopy -d1 -s ./${each_wrf} ./${each_wrf}_compressed 

    # Copy to nextcloud
    rclone copy -vc ./${each_wrf} climatex:Share_Data/WPSDIR/METGRIDDIR/WRFOUT/
    sleep 1
    # Repeat to confirm no change
    rclone copy -vc ./${each_wrf} climatex:Share_Data/WPSDIR/METGRIDDIR/WRFOUT/

    # Check that it's properly on the remote server, and delete if it's there
    rclone ls climatex:Share_Data/WPSDIR/METGRIDDIR/WRFOUT/${each_wrf}_compressed || code=1

    #if [[ $code == 0 ]]; then
    #    rm -v $each_wrf
    #fi

done

rm -fv ./wrfrst*_06*00
rm -fv ./wrfrst*_18*00
rm -fv ./core*

$slack -c $channel -s LOW -t "${rundir} finished wrfout compression" -u $username "$SECONDS s"

touch RCLONE_COMPRESS.OK


