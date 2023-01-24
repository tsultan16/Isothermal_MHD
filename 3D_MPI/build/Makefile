# Compiler
MPIFC = mpifort  # mpiifort 
FC = gfortran  

# The CFLAGS variable sets compile flags:
#  -g                            compile with debug information
#  -Wall                         give verbose compiler warnings
#  -O0                           do not optimize generated code
#  -O2 -ftree-vectorize          apply level 2 optimization with auto-vectorization only
#  -O3                           apply level 3 optimizations
#  -std=gnu99                    use the GNU99 standard language definition
#  -fopt-info-all=opt_report.all save optimization report to file
# -fopenmp                       enable openmp
# -fcheck=all                    enable array bounds checking
# -CB                            intel compiler array bounds checkingg
# -ffpe-trap=invalid             stops program at first NaN (gnu compiler)
# -fpe0                          stops program at first NaN (intel compiler)

CFLAGS = -O3 # -fpe0 #-fopenmp  

OBJS = constants_mod.o grid_data_mod.o isotvd_mod.o vanleer_mod.o io_mod.o  \
       mpi_domain_mod.o init_domain_mod.o boundary_conditions_mod.o turb_driver_mod.o \
       passive_mod.o driver.o 

isotvd_3d : $(OBJS)
	$(MPIFC) $(CFLAGS) -o ../isotvd_3d $(OBJS)
	
constants_mod.o : ../src/constants_mod.f90
	$(MPIFC) $(CFLAGS) -c ../src/constants_mod.f90
		
grid_data_mod.o : ../src/grid_data_mod.f90 constants_mod.o
	$(MPIFC) $(CFLAGS) -c ../src/grid_data_mod.f90	

init_domain_mod.o : ../src/init_domain_mod.f90  constants_mod.o grid_data_mod.o  
	$(MPIFC) $(CFLAGS) -c ../src/init_domain_mod.f90

boundary_conditions_mod.o : ../src/boundary_conditions_mod.f90  constants_mod.o grid_data_mod.o  
	$(MPIFC) $(CFLAGS) -c ../src/boundary_conditions_mod.f90

turb_driver_mod.o : ../src/turb_driver_mod.f90  constants_mod.o grid_data_mod.o mpi_domain_mod.o 
	$(MPIFC) $(CFLAGS) -c ../src/turb_driver_mod.f90
    
isotvd_mod.o : ../src/isotvd_mod.f90  constants_mod.o grid_data_mod.o  
	$(MPIFC) $(CFLAGS) -c ../src/isotvd_mod.f90

passive_mod.o : ../src/passive_mod.f90  constants_mod.o grid_data_mod.o  
	$(MPIFC) $(CFLAGS) -c ../src/passive_mod.f90
    
vanleer_mod.o : ../src/vanleer_mod.f90  constants_mod.o grid_data_mod.o turb_driver_mod.o 
	$(MPIFC) $(CFLAGS) -c ../src/vanleer_mod.f90
    
mpi_domain_mod.o : ../src/mpi_domain_mod.f90  constants_mod.o grid_data_mod.o
	$(MPIFC) $(CFLAGS) -c ../src/mpi_domain_mod.f90

io_mod.o : ../src/io_mod.f90  constants_mod.o grid_data_mod.o mpi_domain_mod.o
	$(MPIFC) $(CFLAGS) -c ../src/io_mod.f90	

    
driver.o : ../src/driver.f90  constants_mod.o grid_data_mod.o io_mod.o mpi_domain_mod.o \
                           isotvd_mod.o vanleer_mod.o turb_driver_mod.o passive_mod.o
	$(MPIFC) $(CFLAGS) -c ../src/driver.f90
	    
    
clean:
	rm  -f *.o *.mod *.txt