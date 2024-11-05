#!/bin/bash -l

set -x

module load ips/18.0.1.163 grib_util/1.0.6 prod_envir/1.0.3 prod_util/1.1.4 NetCDF/4.5.0 HDF5-serial/1.10.1

HOST1=vm-lnx-cpcwork1
HOST2=vm-lnx-cpcwork2

if test $# -eq 1
then
PDY=$1
else
PDY=`./cal_date_y4`
fi

PDYm1=`./cal_date_y4 -1 ${PDY}`
PDYm3=`./cal_date_y4 -10 ${PDY}`

target=`echo $PDY | cut -c3-8`
old=`echo $PDYm3 | cut -c3-8`
rm -f $target.fcst $old.fcst

#COM=$COMROOT
COM=`compath.py gfs/prod`

#######################################################################
# Setup environment variable XLFRTEOPS to use unit_vars
#######################################################################
export XLFRTEOPTS="unit_vars = YES"
export XLFUNIT_52="${target}.fcst"

rm -f fort.52 ${target}.fcst

#######################################################################
# Get land-sea mask
#######################################################################
rm -f landmask_1536.dat
$WGRIB2 ${COM}/gfs.${PDY}/00/atmos/gfs.t00z.sfluxgrbf000.grib2 | grep ":LAND:" | $WGRIB2 -i ${COM}/gfs.${PDY}/00/atmos/gfs.t00z.sfluxgrbf000.grib2 -order we:ns -bin landmask_1536.dat

#######################################################################
# Calculate 12 Hour Forecast
#######################################################################
for FHOUR in 012 024 036 048 060 072 084 096 108 120 132 144 156 168 180
do
   ln -sf ${COM}/gfs.${PDY}/00/atmos/gfs.t00z.atmf${FHOUR}.nc ./drfmr.t00z.sf
   export XLFUNIT_11="drfmr.t00z.sf"
   ./fcst_driver_nc.x drfmr.t00z.sf
   rm -f drfmr.t00z.sf
done

rm -f landmask_1536.dat

mv fort.52 "${target}.fcst"

ssh wd53rl@cpcrzdm bin/clean_ftp.sh

filelen=`cat ${target}.fcst | wc -c`
if test $filelen -lt 2310
then
rm -f ${target}.fcst
exit
fi

scp -p ${target}.fcst wd53rl@cpcrzdm:ftp/new_sbaam/.
scp -p ${target}.fcst roger.lin@${HOST1}:new_sbaam/.
scp -p ${target}.fcst roger.lin@${HOST2}:new_sbaam/.

find . -maxdepth 1 -name '*.fcst' -type f -mtime +14
find . -maxdepth 1 -name '*.fcst' -type f -mtime +14 -exec rm -f {} \;

exit
