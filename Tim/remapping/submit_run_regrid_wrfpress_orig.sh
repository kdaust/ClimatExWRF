#!/bin/bash
#SBATCH --time=23:00:00
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

rundir=RUNDIR

$slack -c $channel -s MED -t "${rundir} started wrfpress processing" -u $username "JOBNAME"

cd $rundir

if [[ ! -e WRFPRESS.OK ]]; then
    time ./run_regrid_wrfpress.sh >& WRFPRESS.log
fi

if [[ ! -e WRFPRESS.OK ]]; then
    exit 1
fi

rclone mkdir climatex:Share_Data/WPSDIR/METGRIDDIR/WRFPRESS
SECONDS=0

# WRFPRESS files
for each_wrf in `ls ./wrfpress*00`; do
    code=0
    echo $each_wrf

    # Copy to nextcloud
    rclone copy -vc ./${each_wrf} climatex:Share_Data/WPSDIR/METGRIDDIR/WRFPRESS/ --max-depth=1
    sleep 1
    # Repeat to confirm no change
    rclone copy -vc ./${each_wrf} climatex:Share_Data/WPSDIR/METGRIDDIR/WRFPRESS/ --max-depth=1

    # Check that it's properly on the remote server, and delete if it's there
    rclone ls climatex:Share_Data/WPSDIR/METGRIDDIR/WRFPRESS/${each_wrf} || code=1 

    #if [[ $code == 0 ]]; then
    #    rm -v $each_wrf
    #fi
    nccopy -d1 -s ./${each_wrf} ./${each_wrf}_compressed
    rm -fv ./${each_wrf}
done

# GHT files
for each_wrf in `ls ./GHT*`; do
    code=0
    echo $each_wrf

    # Copy to nextcloud
    rclone copy -vc ./${each_wrf} climatex:Share_Data/WPSDIR/METGRIDDIR/WRFPRESS/ --max-depth=1
    sleep 1
    # Repeat to confirm no change
    rclone copy -vc ./${each_wrf} climatex:Share_Data/WPSDIR/METGRIDDIR/WRFPRESS/ --max-depth=1

    # Check that it's properly on the remote server, and delete if it's there
    rclone ls climatex:Share_Data/WPSDIR/METGRIDDIR/WRFPRESS/${each_wrf} || code=1


done
$slack -c $channel -s LOW -t "${rundir} finished wrfpress processing" -u $username "$SECONDS s"
touch RCLONE_WRFPRESS.OK

