#!/bin/bash
#SBATCH --time=15:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem-per-cpu=3950M
#SBATCH --job-name=PGWJOB
#SBATCH --account=rrg-rstull


ulimit -s unlimited
ulimit -v unlimited
ulimit -m unlimited

user=`whoami`
slack=/home/${user}/climatex/slack.sh
username="$(whoami)@$(hostname -s)"
channel=climatex_runs
target='<!channel>'

$slack -c $channel -s MED -t "PGWJOB STARTED" -u $username "$target PIDIR"

source /home/${user}/WRFV4.3.3/source_wrf

code=0
./pgw_era5.sh YYYY MM SKIP || code=1

if [[ $code == 1 ]]; then
    $slack -c $channel -s HIGH -t "PGWJOB FAILED" -u $username "$target PIDIR"
    exit 1
else
    $slack -c $channel -s LOW -t "PGWJOB DONE" -u $username "$target PIDIR"
    touch PGWJOB.OK
    touch ../PGWJOB.OK

fi
