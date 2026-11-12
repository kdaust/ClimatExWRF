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

$slack -c $channel -s MED -t "${rundir} started wrf gzip" -u $username "JOBNAME"

cd $rundir
    
    metgrid_dir=${each_dir}
    if [ -e REAL.OK ]; then
	for each_wrfinput in `ls wrf*_d0{1,2,3}`; do
	    # Copy wrfinputs, wrfbdy and wrfdda to nextcloud
	    rclone copy -vc ./${each_wrfinput} climatex:Share_Data/WPSDIR/METGRIDDIR/
	    sleep 1
	done
    fi


for each_remap in `ls | grep REMAP`; do
    code=0
    echo $each_remap
    
    tar -cf ZIPPED_${each_remap}.zip ${each_remap}
    rclone copy -vc ZIPPED_${each_remap}.zip climatex:Share_Data/WPSDIR/METGRIDDIR/

    # Check that it's properly on the remote server, and delete if it's there
    rclone ls climatex:Share_Data/WPSDIR/METGRIDDIR/ZIPPED_${each_remap}.zip || code=1

done


$slack -c $channel -s LOW -t "${rundir} finished wrf gzip" -u $username "$SECONDS s"

touch RCLONE_GZIP.OK


