PROGRAM isothermal_tvd_driver

USE constants_mod
USE grid_data_mod
USE Iso_TVD_Solver_mod
USE io_mod
USE mpi_domain_mod
USE boundary_conditions_mod
USE init_domain_mod

IMPLICIT NONE

INTEGER :: tsteps, t_begin, t_end
REAL*8  :: reduced_dt
INTEGER :: reduce_req, parallel_io_req, parallel_file_handle
REAL*8  :: ttot = 0.d0, tsolve = 0.d0, tfile = 0.d0, tcom = 0.d0, trma = 0.d0, tred = 0.d0, tsolve_1, tsolve_2,  tfile_1, tfile_2, t1, t2, tcom_1, tcom_2, tcom_3, tcom_4, &
           tcom_5, tcom_6, tcom_7, tcom_8, tcom_9, tcom_10, tfile_3, tfile_4   
   
   
! start MPI and set up domain decomposition
IF(myrank .EQ. 0) PRINT*,'Initializing MPI...'
CALL initialize_MPI()

! allocate memory for grid data
CALL create_grid_arrays()

! initialize solver
CALL init_solver()

! initialize fluid state
IF(.NOT. read_from_restart_file) THEN
    
    CALL init_state()
    t_begin = 1
    t_end = maxsteps
    
    ! Exchange boundary data with neighbors
    IF(myrank .EQ. 0) PRINT*,'Exchanging boundary for initial state.'
    CALL initiate_halo_exchange()
    CALL exchange_bndry_data()
    CALL end_halo_exchange()
    IF(myrank .EQ. 0) PRINT*,'Boundary exchange completed.'

    CALL boundary_conditions()
    
ELSE
    
    CALL read_restart(t_begin)
    t_end = t_begin + maxsteps - 1
    
END IF

IF(.NOT. parallel_io) THEN
    CALL writetofile_unformatted(t_begin-1)
ELSE
    CALL load_parallel_file_buffer()
    CALL writetofile_unformatted_parallel(t_begin-1, parallel_io_req, parallel_file_handle, .TRUE., .FALSE.)
END IF


! set initial time step size
my_dt = 1.d-15 


CALL MPI_BARRIER(comm3d, ierr)

t1 = MPI_Wtime()

! simulation loop
DO tsteps = t_begin, t_end

    IF(myrank .EQ. 0) PRINT*,'Time step = ',tsteps
    IF(myrank .EQ. 0) PRINT*,'dt =', my_dt
    IF(myrank .EQ. 0) PRINT*,'% Complete = ', 100.0*(tsteps-t_begin)/(1.0*(t_end-t_begin))

    t2 = MPI_Wtime()
    IF(myrank .EQ. 0) WRITE(*,'("Time elapsed = ",i3," hour ",i3," min ",i3, "seconds.")') INT(t2-t1)/3600 , MOD(INT(t2-t1),3600)/60, MOD(INT(t2-t1),60) 
    IF(myrank .EQ. 0) PRINT*,''
    
    ! initiate halo exchange with neighbors (i.e. send signal to neighbors telling them that we are ready to send)
    IF(myrank .EQ. 0 .AND. print_debug) PRINT*,'Initiating halo exchange.'
    tcom_1 = MPI_Wtime()
    CALL initiate_halo_exchange()
    tcom_2 = MPI_Wtime()

    ! evolve state vector
    tsolve_1 = MPI_Wtime()
    CALL solve_3d(my_dt)
    tsolve_2 = MPI_Wtime()
    
    tsolve = tsolve + tsolve_2 - tsolve_1
    
    ! do the PUTs into neighbor's boundary
    tcom_9 = MPI_Wtime()
    CALL exchange_bndry_data()
    tcom_10 = MPI_Wtime()
    
    ! can do other work here instead of waiting for halo exchange to complete...
    tfile_1 = MPI_Wtime()
    IF(parallel_io .AND. MOD(tsteps,tSkip) .EQ. 0) CALL load_parallel_file_buffer()
    tfile_2 = MPI_Wtime()
    
    
    ! wait for halo exchange to complete
    IF(myrank .EQ. 0 .AND. print_debug) PRINT*,'Completing halo exchange.'
    tcom_3 = MPI_Wtime()
    CALL end_halo_exchange()
    tcom_4 = MPI_Wtime()
    
    ! apply boundary conditions (make sure to do this after halo exchange with neighbor's is done)
    CALL boundary_conditions()
 
    ! compute new local time step
    CALL compute_new_time_step(my_dt)    


    ! save state to file 
    IF(MOD(tsteps,tSkip) .EQ. 0) THEN
    
        tfile_3 = MPI_Wtime()
        IF(print_debug) PRINT*,'Writing to file...'
        
        IF(.NOT. parallel_io) THEN
            CALL writetofile_unformatted(tsteps)
        ELSE
            CALL writetofile_unformatted_parallel(tsteps, parallel_io_req, parallel_file_handle, .FALSE., .FALSE.)
        END IF       
   
        IF(print_debug)PRINT*,'Done writing to file..'
        
        ! also save to restart file if we need to...
        
        tfile_4 = MPI_Wtime()
        
    END IF    
    tfile = tfile + tfile_2 - tfile_1 + tfile_4 - tfile_3

    
    ! post a non-blocking mpi all reduce to communicate my time step
    ! value with all other ranks    
    IF(myrank .EQ. 0 .AND. print_debug) PRINT*,'Initiating mpi_all_reduce.'
    tcom_5 = MPI_Wtime()
    CALL share_dt() 
    tcom_6 = MPI_Wtime()
    
      
    ! maybe do some other work here instead of just waiting for allreduce...
    
    
    ! wait for MPI reduction to complete and get the reduced time step value
    IF(myrank .EQ. 0 .AND. print_debug) PRINT*,'Waiting for all_reduce to complete...'
    tcom_7 = MPI_Wtime()
    CALL get_new_dt()
    tcom_8 = MPI_Wtime()
    
    trma = trma + tcom_2-tcom_1 + tcom_4-tcom_3 + tcom_10-tcom_9
    tred = tred + tcom_6-tcom_5 + tcom_8-tcom_7

END DO

! wait for final parallel file write to complete
IF(parallel_io ) CALL writetofile_unformatted_parallel(tsteps, parallel_io_req, parallel_file_handle, .FALSE., .TRUE.)


t2 = MPI_Wtime()

tcom = trma + tred
ttot = t2-t1

CALL MPI_BARRIER(comm3d, ierr)

! write to restart file
CALL write_restart(tsteps)

CALL MPI_BARRIER(comm3d, ierr)

CALL get_timing_stats()

CALL MPI_BARRIER(comm3d, ierr)


! shut off MPI
CALL shut_off_MPI()

! ********************************
CALL destroy_solver()
CALL destroy_grid_arrays()



!PRINT*,''
!WRITE(*, FMT = '("Myrank = ",i4,", Total execution time (s) = ",f10.5,", Solve time fraction = ",f5.4,", Total MPI comm. time fraction = ",f5.4, &
!      ", RMA time fraction = ",f8.6,", MPI_ALLREDUCE time fraction = ",f8.6,",I/O time fraction =",f8.7)') myrank, ttot, tsolve/ttot, tcom/ttot,&
!        trma/tcom, tred/tcom, tfile/ttot
!PRINT*,''



IF(myrank .EQ. 0) PRINT*,'Done!'

CONTAINS


SUBROUTINE share_dt()

   ! initiate non-blocking MPI all_reduce
   CALL start_allreduce(my_dt, reduced_dt, reduce_req)

END SUBROUTINE share_dt



SUBROUTINE get_new_dt() 

   ! wait for non-blocking MPI all_reduce to complete
   CALL end_allreduce(reduce_req) 
  
   ! replace current time step with the global minimum 
   my_dt = reduced_dt       

END SUBROUTINE get_new_dt


SUBROUTINE get_timing_stats()

    REAL*8 :: all_tot, all_tsolve, all_tcom, all_tfile, t_per_step 
    INTEGER :: ierr

    IF(myrank .EQ. 0) PRINT*,'Getting more timing stats...'

    CALL MPI_ALLREDUCE(ttot, all_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
    CALL MPI_ALLREDUCE(tsolve, all_tsolve, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
    CALL MPI_ALLREDUCE(tcom, all_tcom, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
    CALL MPI_ALLREDUCE(tfile, all_tfile, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)

    IF(myrank .EQ. 0) THEN
    
        PRINT*,''
        PRINT*,'Average per rank:'
    
        WRITE(*, FMT = '("Total execution time (s) = ",f10.5,", Solve time fraction = ",f8.7,", Total MPI comm. time fraction = ",f8.7, &
              ", Total I/O time fraction = ",f8.7)') all_tot/np, all_tsolve/all_tot, all_tcom/all_tot, all_tfile/all_tot
    
        PRINT*,''
        
        OPEN(UNIT = 10, FILE = 'Output/Run_Stats/timing_stats.txt')
        
        WRITE(10, FMT= '("Total execution time (s) = ",f10.2)') all_tot/np
        WRITE(10, FMT= '("Solve Time Fraction = ",f15.14)') all_tsolve/all_tot
        WRITE(10, FMT= '("Halo Exchange Time Fraction = ",f15.14)') all_tcom/all_tot
        WRITE(10, FMT= '("I/O Time Fraction = ",f15.14)') all_tfile/all_tot
               
        CLOSE(UNIT = 10)
        
    END IF

END SUBROUTINE get_timing_stats



END PROGRAM isothermal_tvd_driver