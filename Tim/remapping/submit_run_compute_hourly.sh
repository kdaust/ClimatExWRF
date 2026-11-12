#!/bin/bash
#SBATCH --time=04:00:00
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

domain=$1

source /home/${user}/WRFV4.3.3/source_wrf

module load cdo/2.0.5

rundir=RUNDIR

cd $rundir 

code=0

SECONDS=0
python3 /home/${user}/climatex/remapping/compute_hourly.py $domain || code=1
duration=$SECONDS

if [ $code == 0 ]; then
    $slack -c $channel -s LOW -t "JOBNAME has finished run_remap_wrfhourly" -u $username "$rundir"
    touch WRFHOURLY.${domain}.OK
    exit 0
else
    $slack -c $channel -s HIGH -t "JOBNAME has failed to complete run_remap_wrfhourly" -u $username "$rundir"
    exit 1
fi



