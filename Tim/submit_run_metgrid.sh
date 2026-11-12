#!/bin/bash
#SBATCH --time=11:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem-per-cpu=3950M
#SBATCH --job-name=METJOB
#SBATCH --account=def-rstull


ulimit -s unlimited
ulimit -v unlimited
ulimit -m unlimited

user=`whoami`
slack=/home/${user}/climatex/slack.sh
username="$(whoami)@$(hostname -s)"
channel=climatex_runs
target='<!channel>'

$slack -c $channel -s MED -t "METJOB STARTED" -u $username "$target METDIR"

source /home/${user}/WRFV4.3.3/source_wrf

cd METDIR

SECONDS=0
time ./metgrid.exe 
duration=$SECONDS
code=$?

if [ $code == 0 ]; then
    $slack -c $channel -s LOW -t "METJOB finished in ${duration} s" -u $username "$target METDIR"
    touch METGRID.OK
    echo "METGRID SUCCESS"
else
    $slack -c $channel -s HIGH -t "METJOB FAILED" -u $username "$target METDIR"
    echo "METGRID FAILED"
fi





