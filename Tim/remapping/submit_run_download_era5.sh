#!/bin/bash


ulimit -s unlimited
ulimit -v unlimited
ulimit -m unlimited

download_directory=DOWNLOAD_DIRECTORY
download_year=DOWNLOAD_YEAR
download_month=DOWNLOAD_MONTH
download_land=DOWNLOAD_LAND

cd ${download_directory}

user=`whoami`
slack=/home/${user}/climatex/slack.sh
username="$(whoami)@$(hostname -s)"
channel=climatex_runs
target='<!channel>'

$slack -c $channel -s MED -t "DJOB STARTED" -u $username "$target $download_directory"

source /home/${user}/WRFV4.3.3/source_wrf

SECONDS=0
time python3 /home/${user}/climatex/run_download_era5.py $download_directory $download_year $download_month $download_land
duration=$SECONDS
code=$?

if [ $code == 0 ]; then
    $slack -c $channel -s LOW -t "DJOB finished in ${duration} s" -u $username "$target $download_directory"
    touch IC.OK
    echo "Downloading ERA5 SUCCESS"
    mv -v ERA5* ../
else
    $slack -c $channel -s HIGH -t "DJOB FAILED" -u $username "$target $download_directory"
    echo "Downloading ERA5 FAILED"
fi





