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

$slack -c $channel -s MED -t "${rundir} started wrfuvic processing" -u $username "JOBNAME"

cd $rundir

if [[ ! -e WRFUVIC.OK ]]; then
    time ./run_regrid_wrfuvic.sh >& WRFUVIC.log
fi

if [[ ! -e WRFUVIC.OK ]]; then
    exit 1
fi

rclone mkdir climatex:Share_Data/WPSDIR/METGRIDDIR/WRFUVIC
SECONDS=0

# Uvic files
for each_wrf in `ls ./wrfuvic*`; do
    code=0
    echo $each_wrf

    # Copy to nextcloud
    rclone copy -vc ./${each_wrf} climatex:Share_Data/WPSDIR/METGRIDDIR/WRFUVIC/ --max-depth=1
    sleep 1
    # Repeat to confirm no change
    rclone copy -vc ./${each_wrf} climatex:Share_Data/WPSDIR/METGRIDDIR/WRFUVIC/ --max-depth=1

    # Check that it's properly on the remote server, and delete if it's there
    rclone ls climatex:Share_Data/WPSDIR/METGRIDDIR/WRFUVIC/${each_wrf} || code=1

    #if [[ $code == 0 ]]; then
    #    rm -v $each_wrf
    #fi

done

# Xtrm files
for each_wrf in `ls ./wrfxtrm*`; do
    code=0
    echo $each_wrf

    # Copy to nextcloud
    rclone copy -vc ./${each_wrf} climatex:Share_Data/WPSDIR/METGRIDDIR/WRFUVIC/ --max-depth=1
    sleep 1
    # Repeat to confirm no change
    rclone copy -vc ./${each_wrf} climatex:Share_Data/WPSDIR/METGRIDDIR/WRFUVIC/ --max-depth=1

    # Check that it's properly on the remote server, and delete if it's there
    rclone ls climatex:Share_Data/WPSDIR/METGRIDDIR/WRFUVIC/${each_wrf} || code=1

    #if [[ $code == 0 ]]; then
    #    rm -v $each_wrf
    #fi

done
$slack -c $channel -s LOW -t "${rundir} finished wrfuvic processing" -u $username "$SECONDS s"
touch RCLONE_WRFUVIC.OK

