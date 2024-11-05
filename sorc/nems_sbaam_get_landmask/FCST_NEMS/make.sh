#!/bin/sh
#set -x

mac=$(hostname -f)

case $mac in

#---------------------------------------------------------------------------------
#  WCOSS CRAY.
#---------------------------------------------------------------------------------
llogin? | slogin?)

  module purge
  module load modules/3.2.6.7
  module load PrgEnv-intel/5.2.56
  module rm intel
  module load intel/16.3.210
  module load cray-mpich/7.2.0
  module load craype-haswell
  module load cray-netcdf
  module load w3nco-intel/2.0.6

  export NETCDFPATH="/opt/cray/netcdf/4.3.2/INTEL/140"
  export NETCDF_LIB="-L${NETCDFPATH}/lib -lnetcdf -lnetcdff"
  export NETCDF_INC="-I${NETCDFPATH}/include"
  export FCOMP=ftn
  export FFLAGS="-O3"

  make clean
  make
  rc=$?  ;;

#---------------------------------------------------------------------------------
# THEIA.
#---------------------------------------------------------------------------------

tfe??)

  source /apps/lmod/lmod/init/sh
  module purge

  module load intel/15.1.133
  module load impi/5.1.1.109
  module load netcdf/4.3.0

  export NETCDFPATH="/apps/netcdf/4.3.0-intel"
  export NETCDF_LIB="-L${NETCDFPATH}/lib -lnetcdf -lnetcdff"
  export NETCDF_INC="-I${NETCDFPATH}/include"
  export FCOMP=ifort
  export FFLAGS="-O3"

  make clean
  make
  rc=$?  ;;

#---------------------------------------------------------------------------------
# HERA.
#---------------------------------------------------------------------------------

hfe??)

  source /apps/lmod/lmod/init/sh
  module purge

  module load intel/18.0.5.274 
  module load impi/2018.0.4 
  module load netcdf/4.7.0 
  module load hdf5/1.10.5

  export NETCDFPATH="/apps/netcdf/4.7.0/intel/18.0.5.274"
  export NETCDF_LIB="-L${NETCDFPATH}/lib -lnetcdf -lnetcdff "
  export NETCDF_INC="-I${NETCDFPATH}/include"
  export FCOMP=ifort
  export FFLAGS="-O3"

  make clean
  make
  rc=$?  ;;

#---------------------------------------------------------------------------------
# ON DELL.
#---------------------------------------------------------------------------------

m????.ncep.noaa.gov | v????.ncep.noaa.gov )

  module purge
  module use /usrx/local/dev/modulefiles

  module load ips/18.0.1.163
  module load impi/18.0.1
  module load NetCDF/4.5.0
  module load HDF5-serial/1.10.1
  module load sp/2.0.2

# export NETCDFPATH="/usrx/local/prod/packages/ips/18.0.1/netcdf/4.5.0"
# export NETCDF_LIB="-L${NETCDFPATH}/lib -lnetcdf -lnetcdff"
# export NETCDF_INC="-I${NETCDFPATH}/include"
# export FCOMP=ifort
# export FFLAGS="-O3"

# make clean
  make
  rc=$? ;;

esac

exit
