#!/bin/bash -l

cd /gpfs/dell2/cpc/noscrub/${USER}/nems_sbaam_get_landmask

HOST=`hostname`
HOST_1=`echo $HOST | cut -c1-1`
PROD_1=`cat /etc/prod | cut -c1-1`

if test ${HOST_1}"X" != ${PROD_1}"X"
then
   echo "DEV machine, no processing delay"
else
   echo "PROD machine, delay 10 minutes so it finishes after DEV run"
   sleep 600
fi

./anal_energy_nc.sh 1> anal_energy_nc.stdout 2> anal_energy_nc.stderr
./fcst_energy_nc.sh 1> fcst_energy_nc.stdout 2> fcst_energy_nc.stderr

exit
