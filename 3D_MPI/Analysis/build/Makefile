# Compiler
MPIFC = mpifort   
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

CFLAGS = -O3   

OBJS = constants_mod.o fft_mod.o readfile_mod.o analysis_main.o

do_analysis : $(OBJS)
	$(MPIFC) $(CFLAGS) -o ../do_analysis $(OBJS)
	
constants_mod.o : ../src/constants_mod.f90
	$(MPIFC) $(CFLAGS) -c ../src/constants_mod.f90
		
fft_mod.o : ../src/fft_mod.f90 constants_mod.o
	$(MPIFC) $(CFLAGS) -c ../src/fft_mod.f90	

readfile_mod.o : ../src/readfile_mod.f90  constants_mod.o  
	$(MPIFC) $(CFLAGS) -c ../src/readfile_mod.f90


analysis_main.o : ../src/analysis_main.f90  constants_mod.o fft_mod.o readfile_mod.o 
	$(MPIFC) $(CFLAGS) -c ../src/analysis_main.f90
    
clean:
	rm  -f *.o *.mod 