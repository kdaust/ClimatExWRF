#!/bin/bash
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem-per-cpu=3950M
#SBATCH --job-name=JOBNAME
#SBATCH --account=rrg-rstull

ulimit -s unlimited
ulimit -v unlimited
ulimit -m unlimited

user=`whoami`
slack=/home/${user}/climatex/slack.sh
username="$(whoami)@$(hostname -s)"
channel=climatex_runs
target='<!channel>'

$slack -c $channel -s MED -t "JOBNAME STARTED" -u $username "$target RUNDIR"

#source /home/${user}/WRFV4.3.3/source_wrf

rundir=RUNDIR

cd $rundir 

code=0

SECONDS=0
time mpirun ${rundir}/real.exe || code=1
duration=$SECONDS

if [ $code == 0 ]; then
    touch ${duration}.REAL.OK
    $slack -c $channel -s LOW -t "JOBNAME has finished real.exe" -u $username "$rundir"
    touch REAL.OK
    rclone copy -vc . climatex:Share_Data/WPSDIR/METGRIDDIR/ --include="wrf*_d0{1,2,3}"
    exit 0
else
    $slack -c $channel -s HIGH -t "JOBNAME has failed to complete real.exe" -u $username "$rundir"
    exit 1
fi



