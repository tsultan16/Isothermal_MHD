# Compiler
MPIFC = mpifort
FC = gfortran  


# The CFLAGS variable sets compile flags for gfortran:
#  -g                            compile with debug information
#  -Wall                         give verbose compiler warnings
#  -O0                           do not optimize generated code
#  -O2 -ftree-vectorize          enable auto-vectorization only
#  -O3                           apply level 3 optimizations
#  -std=gnu99                    use the GNU99 standard language definition
#  -fopt-info-all=opt_report.all save optimization report to file
# -fopenmp                       enable openmp
# -fcheck=all                    enable array bounds checking

CFLAGS = -fcheck=all -ffpe-trap=invalid #-O2 -ftree-vectorize -fopenmp  

OBJS = constants_mod.o grid_data_mod.o isotvd_mod.o driver.o

isotvd_1d : $(OBJS)
	$(FC) $(CFLAGS) -o isotvd_1d $(OBJS)
	
constants_mod.o : src/constants_mod.f90
	$(FC) $(CFLAGS) -c src/constants_mod.f90
		
grid_data_mod.o : src/grid_data_mod.f90 constants_mod.o
	$(FC) $(CFLAGS) -c src/grid_data_mod.f90	

isotvd_mod.o : src/isotvd_mod.f90  constants_mod.o grid_data_mod.o  
	$(FC) $(CFLAGS) -c src/isotvd_mod.f90
	
driver.o : src/driver.f90  constants_mod.o grid_data_mod.o isotvd_mod.o  
	$(FC) $(CFLAGS) -c src/driver.f90
	
clean:
	rm  -f *.o *.mod *.txt