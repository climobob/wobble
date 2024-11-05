#!/bin/bash -l

cd /gpfs/dell2/cpc/noscrub/${USER}/nems_sbaam_get_landmask

if test $# -eq 1
then
TARGET=$1
else
TARGET=`./cal_date_y4`
fi

./anal_energy_nc.sh $TARGET
./fcst_energy_nc.sh $TARGET

exit
