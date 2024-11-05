#!/bin/bash -l

set -x

module load ips/18.0.1.163 grib_util/1.0.6 prod_envir/1.0.3 NetCDF/4.5.0 HDF5-serial/1.10.1

HOST1=vm-lnx-cpcwork1
HOST2=vm-lnx-cpcwork2

if test $# -eq 1
then
PDY=$1
else
PDY=`./cal_date_y4`
fi

PDYm1=`./cal_date_y4 -1 ${PDY}`
PDYm4=`./cal_date_y4 -11 ${PDY}`

target=`echo $PDYm1 | cut -c3-8`
old=`echo $PDYm4 | cut -c3-8`
rm -f $target.anl $old.anl

#COM=$COMROOT
COM=`compath.py gfs/para`

#######################################################################
# Setup environment variable XLFRTEOPS to use unit_vars
#######################################################################
export XLFRTEOPTS="unit_vars = YES"
export XLFUNIT_51="${target}.anl"

rm -f fort.51 ${target}.anl

#######################################################################
# Get land-sea mask
#######################################################################
rm -f landmask_1536.dat
$WGRIB2 ${COM}/gdas.${PDYm1}/12/atmos/gdas.t12z.sfluxgrbf000.grib2 | grep ":LAND:" | $WGRIB2 -i ${COM}/gdas.${PDYm1}/12/atmos/gdas.t12z.sfluxgrbf000.grib2 -order we:ns -bin landmask_1536.dat

#######################################################################
# Calculate T00Z - T18Z for day 1
#######################################################################
for CYCLE in 00 06 12 18
do
ln -sf ${COM}/gdas.${PDYm1}/${CYCLE}/atmos/gdas.t${CYCLE}z.atmanl.nc ./gdas1.sanl
export XLFUNIT_11="gdas1.sanl"
./anal_driver_nc.x gdas1.sanl
rm -f gdas1.sanl
done

#######################################################################
# Calculate T00Z for day 2
#######################################################################
ln -sf ${COM}/gdas.${PDY}/00/atmos/gdas.t00z.atmanl.nc ./gdas1.t00z.sanl
export XLFUNIT_11="gdas1.t00z.sanl"
./anal_driver_nc.x gdas1.t00z.sanl
rm -f gdas1.t00z.sanl

rm -f landmask_1536.dat

mv fort.51 "${target}.anl"

exit

