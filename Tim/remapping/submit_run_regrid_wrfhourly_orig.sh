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

$slack -c $channel -s MED -t "${rundir} started wrfhourly processing" -u $username "JOBNAME"

cd $rundir

if [[ ! -e WRFHOURLY.OK ]]; then
    time ./run_regrid_wrfhourly.sh >& WRFHOURLY.log
fi

if [[ ! -e WRFHOURLY.OK ]]; then
    exit 1
fi

rclone mkdir climatex:Share_Data/WPSDIR/METGRIDDIR/WRFHOURLY
SECONDS=0

for each_wrf in `ls ./wrfhourly*`; do
    code=0
    echo $each_wrf

    # Copy to nextcloud
    rclone copy -vc ./${each_wrf} climatex:Share_Data/WPSDIR/METGRIDDIR/WRFHOURLY/ --max-depth=1
    sleep 1
    # Repeat to confirm no change
    rclone copy -vc ./${each_wrf} climatex:Share_Data/WPSDIR/METGRIDDIR/WRFHOURLY/ --max-depth=1

    # Check that it's properly on the remote server, and delete if it's there
    rclone ls climatex:Share_Data/WPSDIR/METGRIDDIR/WRFHOURLY/${each_wrf} || code=1

    #if [[ $code == 0 ]]; then
    #    rm -v $each_wrf
    #fi

done
$slack -c $channel -s LOW -t "${rundir} finished wrfhourly processing" -u $username "$SECONDS s"
touch RCLONE_WRFHOURLY.OK

