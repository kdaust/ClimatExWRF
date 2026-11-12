#!/bin/bash -l
set -x

if [ $# -lt 6 ]; then
    # TODO: print usage
    echo "Usage: run_ungrib SSYEAR SSMONTH SSDAY EEYEAR EEMONTH EEDAY"
    exit 1
fi


SSYEAR=$1
SSMONTH=$2
SSDAY=$3
EEYEAR=$4
EEMONTH=$5
EEDAY=$6

RUNNAME=$7

user=`whoami`

endyearp1=$(( EEYEAR + 1 ))

ROOT_DIR=/scratch/${user}/CLIMATEX

if [ ! -z "$RUNNAME" ]; then
    ungrib_dir=${ROOT_DIR}/WPS_${SSYEAR}${SSMONTH}_${EEYEAR}${EEMONTH}_${user}_${RUNNAME}
else
    ungrib_dir=${ROOT_DIR}/WPS_${SSYEAR}${SSMONTH}_${EEYEAR}${EEMONTH}_${user}
fi

IBCS_DIR=${ungrib_dir}/ibcs/

echo ""
echo "==================="
echo ""
echo "Making directory ${ungrib_dir}"
echo ""
echo "==================="
echo ""

mkdir ${ungrib_dir}

cd ${ungrib_dir} || exit 1

cp /home/${user}/WRFV4.3.3/WPS-4.3.1/link_grib.csh .
ln -s /home/${user}/WRFV4.3.3/WPS-4.3.1/ungrib.exe .

cp /home/${user}/WRFV4.3.3/WPS-4.3.1/ungrib/Variable_Tables/Vtable.ERA-interim.pl ./Vtable
cp /home/${user}/WRFV4.3.3/WPS-4.3.1/ungrib/Variable_Tables/Vtable.ERA-interim.snow ./Vtable.SNOW

cp /home/${user}/climatex/namelist.ungrib ./namelist.ungrib

months=( 01 02 03 04 05 06 07 08 09 10 11 12 )
startmonth_ind=$(( 10#${SSMONTH} - 1 ))

list_string="${SSMONTH}"
startmonth_ind=$(( startmonth_ind + 1 ))
while [ $startmonth_ind -le 11 ]; do
    list_string="${list_string},${months[$startmonth_ind]}"
    startmonth_ind=$(( startmonth_ind + 1 ))
done


# First process ERA5
rm -f ./GRIBFILE*
rm -f ./PFILE*
./link_grib.csh ${IBCS_DIR}/ERA5_{single,pressure}_${SSYEAR}_{${list_string}}* ${IBCS_DIR}/ERA5_{single,pressure}_${EEYEAR}* ${IBCS_DIR}/ERA5_{single,pressure}_${endyearp1}*_01*

cp namelist.ungrib namelist.wps

sed -i -e "s/SSYEAR/$SSYEAR/g" namelist.wps
sed -i -e "s/SSMONTH/$SSMONTH/g" namelist.wps
sed -i -e "s/SSDAY/$SSDAY/g" namelist.wps
sed -i -e "s/EEYEAR/$EEYEAR/g" namelist.wps
sed -i -e "s/EEMONTH/$EEMONTH/g" namelist.wps
sed -i -e "s/EEDAY/$EEDAY/g" namelist.wps

sed -i -e "s/UNGRIB_FORMAT/FILE/g" namelist.wps

./ungrib.exe || exit 1

# Compute TAVGSFC
ln -s /home/${user}/WRFV4.3.3/WPS-4.3.1/util/avg_tsfc.exe .
./avg_tsfc.exe || exit 1

# Then link in LAND
rm -f ./GRIBFILE*
./link_grib.csh ${IBCS_DIR}/ERA5_land_${SSYEAR}_${SSMONTH}* 

ssdayp1=`printf %02d $(( 10#${SSDAY} + 1 ))`

cp namelist.ungrib namelist.wps
sed -i -e "s/SSYEAR/$SSYEAR/g" namelist.wps
sed -i -e "s/SSMONTH/$SSMONTH/g" namelist.wps
sed -i -e "s/SSDAY/$SSDAY/g" namelist.wps
sed -i -e "s/EEYEAR/$SSYEAR/g" namelist.wps
sed -i -e "s/EEMONTH/$SSMONTH/g" namelist.wps
sed -i -e "s/EEDAY/$ssdayp1/g" namelist.wps

sed -i -e "s/UNGRIB_FORMAT/LAND/g" namelist.wps

./ungrib.exe || exit 1

# Handle SNOW
rm -f ./GRIBFILE*
./link_grib.csh ${IBCS_DIR}/ERA5_single_${SSYEAR}_${SSMONTH}* 

cp namelist.ungrib namelist.wps
sed -i -e "s/SSYEAR/$SSYEAR/g" namelist.wps
sed -i -e "s/SSMONTH/$SSMONTH/g" namelist.wps
sed -i -e "s/SSDAY/$SSDAY/g" namelist.wps
sed -i -e "s/EEYEAR/$SSYEAR/g" namelist.wps
sed -i -e "s/EEMONTH/$SSMONTH/g" namelist.wps
sed -i -e "s/EEDAY/$ssdayp1/g" namelist.wps

sed -i -e "s/UNGRIB_FORMAT/SNOW/g" namelist.wps
cp Vtable.SNOW Vtable

./ungrib.exe || exit 1

touch UNGRIB.OK 






