SHELL=/bin/sh
FCOMP=ifort
FFLAGS="-O3"

CMD=fcst_driver_nc.x

OBJS=fcst_driver_nc.o fcst_energy.o date_offset.o

ALL: $(CMD)

fcst_driver_nc.x: $(OBJS)
	 $(FCOMP) $(FFLAGS)  -o $@ $(OBJS) -I$(NETCDF_INC) -L$(NETCDF_LIB) -lnetcdf -lnetcdff $(SP_LIB4)

fcst_driver_nc.o:  fcst_driver_nc.f90
	 $(FCOMP) $(FFLAGS) -I$(NETCDF_INC) -L$(NETCDF_LIB) -lnetcdf -lnetcdff $(SP_LIB4) -c fcst_driver_nc.f90

fcst_energy.o: fcst_energy.f
	$(FCOMP) -FR -c fcst_energy.f

date_offset.o: date_offset.f
	$(FCOMP) -FR -c date_offset.f

clean:
	 rm -f *.o *.mod ${CMD} 
