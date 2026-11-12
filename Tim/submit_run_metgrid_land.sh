#!/bin/bash
#SBATCH --time=03:45:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=3950M
#SBATCH --constraint="[skylake|cascade]"
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
    echo "METGRID SUCCESS"
else
    $slack -c $channel -s HIGH -t "METJOB FAILED" -u $username "$target METDIR"
    echo "METGRID FAILED"
    exit 1
fi


module load python/3.10.2
module load mpi4py/3.1.3
module load  geos/3.9.1
module load  gdal/3.0.4
module load proj/9.0.1

source /home/${user}/cartopy/bin/activate

year=SSYEAR
month=SSMONTH

for each_met in `ls met_em*${year}-${month}-01*`; do
    cp -v ./${each_met} ./era5_${each_met}
    python3 /home/${user}/climatex/merge_metgrid.py ./${each_met} ./land_${each_met} || exit 1
done

touch METGRID.OK
touch METGRID_LAND.OK





