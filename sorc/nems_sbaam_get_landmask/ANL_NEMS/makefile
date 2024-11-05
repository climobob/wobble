SHELL=/bin/sh
FCOMP=ifort
FFLAGS="-O3"

CMD=anal_driver_nc.x

OBJS=anal_driver_nc.o anal_energy.o

ALL: $(CMD)

anal_driver_nc.x: $(OBJS)
	 $(FCOMP) $(FFLAGS)  -o $@ $(OBJS) -I$(NETCDF_INC) -L$(NETCDF_LIB) -lnetcdf -lnetcdff $(SP_LIB4)

anal_driver_nc.o:  anal_driver_nc.f90
	 $(FCOMP) $(FFLAGS) -I$(NETCDF_INC) -L$(NETCDF_LIB) -lnetcdf -lnetcdff $(SP_LIB4) -c anal_driver_nc.f90

anal_energy.o: anal_energy.f
	$(FCOMP) -FR -c anal_energy.f

clean:
	 rm -f *.o *.mod ${CMD} 
