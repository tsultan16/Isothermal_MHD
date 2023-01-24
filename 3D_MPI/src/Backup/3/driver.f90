PROGRAM isothermal_tvd_driver

USE constants_mod
USE grid_data_mod
USE Iso_TVD_Solver_mod
USE VanLeer_Solver_mod
USE io_mod
USE mpi_domain_mod
USE boundary_conditions_mod
USE init_domain_mod
USE passive_solver_mod
USE turb_driver_mod

IMPLICIT NONE


REAL*8  :: reduced_dt, t_sim
INTEGER :: reduce_req, parallel_io_req, parallel_file_handle
REAL*8 :: my_divB, t_old, t_dump, t_restart, t_turb   
LOGICAL :: write_status, drive_next
   
! do initializations for grid, solver, mpi and i/o to get everything started   
CALL main_init()
    
CALL MPI_BARRIER(comm3d, ierr)

t1 = MPI_Wtime()

! run the simulation loop
DO WHILE((t_sim .LT. t_end) .AND. (tsteps .LT. maxsteps))

    IF(myrank .EQ. 0) PRINT*,'Time step = ',tsteps
    IF(myrank .EQ. 0) PRINT*,'dt =', my_dt
    IF(myrank .EQ. 0) PRINT*,'Simulation Clock Time (code units) =', t_sim
    IF(myrank .EQ. 0) PRINT*,'% Complete = ', 100.0 * t_sim / t_end

    t2 = MPI_Wtime()
    IF(myrank .EQ. 0) WRITE(*,'("Time elapsed = ",i3," hour ",i3," min ",i3, "seconds.")') INT(t2-t1)/3600 , MOD(INT(t2-t1),3600)/60, MOD(INT(t2-t1),60) 
    IF(myrank .EQ. 0) PRINT*,''
    
    t_old = t_sim
    
    ! evolve MHD state vector and passive variables
    tsolve_1 = MPI_Wtime()
    IF(npass .GT. 0) CALL copy_velocity(.FALSE.)
    
    IF(drive_turbulence .AND. (t_turb .GT. t_old) .AND. (t_turb .LE. t_old+my_dt)) THEN
        drive_next = .TRUE.
        ! compute the next turbulence driving time
        t_turb = MAX(t_turb + turb_frequency, t_sim + 1d-2*my_dt)
     END IF   
        
    CALL mhd_update()
    
    IF(npass .GT. 0) THEN
        CALL copy_velocity(.TRUE.)
        IF(print_debug) PRINT*,'Updating passive variables'
        CALL passive_solve_3d(my_dt)
        IF(print_debug) PRINT*,'Done updating passive variables'
    END IF        
    tsolve_2 = MPI_Wtime()    
    tsolve = tsolve + tsolve_2 - tsolve_1
       
    ! advance the simulation clock time
    t_sim = t_sim + my_dt
    tsteps = tsteps + 1
    
    ! do halo ecxhange with neighbor's boundary neighbor's boundary
    tcom_1 = MPI_Wtime()
    IF(myrank .EQ. 0) PRINT*,'Doing halo exchange..'
    CALL start_halo_exchange()
    tcom_2 = MPI_Wtime()
    
    !*****************************************************************************************************    
    ! Halo exchange (via MPI RMA) is taking place in the background. Meanwhile, we can do other work here. 
    !*****************************************************************************************************
    !*****************************************************************************************************
    
    ! generate random force for next time step
    tsolve_1 = MPI_Wtime()
    IF(drive_turbulence .AND. drive_next) THEN
        CALL generate_rand_force()
        drive_next = .FALSE.    
    END IF    
    tsolve_2 = MPI_Wtime()
    tsolve = tsolve + tsolve_2 - tsolve_1
    
    ! save state to file 
    tfile_3 = MPI_Wtime()
    IF((t_dump .GT. t_old) .AND. (t_dump .LE. t_sim)) THEN
        CALL file_dump(1)
    END IF
    
    IF((t_restart .GT. t_old) .AND. (t_restart .LE. t_sim)) THEN
        CALL file_dump(2)
    END IF
    tfile_4 = MPI_Wtime() 
    
    !*****************************************************************************************************
    !*****************************************************************************************************
    tcom_3 = MPI_Wtime()
    CALL end_halo_exchange()
    IF(myrank .EQ. 0) PRINT*,'Halo exchange completed.'
    tcom_4 = MPI_Wtime()

    ! apply boundary conditions (make sure to do this after halo exchange with neighbor's is done)
    CALL boundary_conditions()
       
    ! compute new local time step
    tsolve_1 = MPI_Wtime()
    CALL compute_new_time_step(my_dt)    
    tsolve_2 = MPI_Wtime()
    tsolve = tsolve + tsolve_2 - tsolve_1
    
        
    ! post a non-blocking mpi all reduce to communicate my time step
    ! value with all other ranks    
    IF(myrank .EQ. 0 .AND. print_debug) PRINT*,'Initiating mpi_all_reduce.'
    tcom_3 = MPI_Wtime()
    CALL share_dt() 
    tcom_4 = MPI_Wtime()
    
    !*****************************************************************************************************
    ! maybe do some other work here instead of just waiting for allreduce to complete...
    !*****************************************************************************************************
    !*****************************************************************************************************
    
    ! compute divB
    CALL compute_divB()   
    
    CALL compute_KE()
    
    !*****************************************************************************************************
    !*****************************************************************************************************
    ! wait for MPI reduction to complete and get the reduced time step value
    IF(myrank .EQ. 0 .AND. print_debug) PRINT*,'Waiting for all_reduce to complete...'
    tcom_5 = MPI_Wtime()
    CALL get_new_dt()
    tcom_6 = MPI_Wtime()    


    tfile = tfile + tfile_4 - tfile_3    
    trma = trma + tcom_2-tcom_1 + tcom_4 - tcom_3
    tred = tred + tcom_4-tcom_3 + tcom_6-tcom_5

END DO

IF(t_sim .LT. t_end .AND. myrank .EQ. 0) PRINT*,'Exceeded maximum number of time steps...exiting simulation loop.'


! wait for final parallel file write to complete
tfile_5 = MPI_Wtime()
IF(parallel_io) CALL writetofile_unformatted_parallel(tsteps, parallel_io_req, parallel_file_handle, .FALSE., .TRUE.,write_status)
tfile_6 = MPI_Wtime()
tfile = tfile + tfile_6 - tfile_5

t2 = MPI_Wtime()
tcom = trma + tred
ttot = t2-t1 ! tsolve + tcom + tfile


CALL MPI_BARRIER(comm3d, ierr)

! write to restart file
CALL write_restart(t_sim)

CALL MPI_BARRIER(comm3d, ierr)

CALL get_timing_stats()

CALL MPI_BARRIER(comm3d, ierr)


! shut off MPI
CALL shut_off_MPI()

! ********************************
IF(mhd_solver_type .EQ. 1) THEN
    CALL destroy_solver_tvd()
ELSE IF(mhd_solver_type .EQ. 2) THEN
    CALL destroy_solver_vanleer()
END IF

IF(npass .GT. 0) CALL destroy_passive_solver()

IF(drive_turbulence) CALL destroy_turb_driver()


CALL destroy_grid_arrays()



!PRINT*,''
!WRITE(*, FMT = '("Myrank = ",i4,", Total execution time (s) = ",f10.5,", Solve time fraction = ",f5.4,", Total MPI comm. time fraction = ",f5.4, &
!      ", RMA time fraction = ",f8.6,", MPI_ALLREDUCE time fraction = ",f8.6,",I/O time fraction =",f8.7)') myrank, ttot, tsolve/ttot, tcom/ttot,&
!        trma/tcom, tred/tcom, tfile/ttot
!PRINT*,''

CLOSE(UNIT = 99)

IF(myrank .EQ. 0) PRINT*,'Done!'


CONTAINS


! Top-level routine for initialization of every major component 
SUBROUTINE main_init()

    INTEGER :: s1_in, s2_in, s3_in, ierr
    LOGICAL :: write_status
    CHARACTER(LEN=180) :: filename


    ! start MPI and set up domain decomposition
    CALL initialize_MPI()

    ! allocate memory for grid data
    CALL create_grid_arrays()

    ! set initial time step size
    my_dt = 1.d-15 
    t_sim = 0.d0

    ! set cell size
    my_dx = 1.d0 / DBLE(MAX(nx*nranks_x,ny*nranks_y,nz*nranks_z))

    ! open a log file for file dumps
    filename = TRIM(output_filepath)//TRIM('/Run_Stats/dump_log.txt')
    OPEN(UNIT = 99 , FILE = filename)

    ! initialize solver
    IF(mhd_solver_type .EQ. 1) THEN
        CALL init_solver_tvd()
    ELSE IF(mhd_solver_type .EQ. 2) THEN
        CALL init_solver_vanleer()
    END IF

    IF(npass .GT. 0) CALL init_passive_solver()

    IF(myrank .EQ. 0) THEN
    
        IF(riemann_solver_type .EQ. 1) THEN
            PRINT*,'Riemann solver type = ROE'
        ELSE IF(riemann_solver_type .EQ. 2) THEN
            PRINT*,'Riemann solver type = HLLD'
        ELSE IF(riemann_solver_type .EQ. 3) THEN
            PRINT*,'Riemann solver type = HLLE'
        END IF
       
    END IF

    CALL MPI_BARRIER(comm3d, ierr)

    ! initialize fluid state
    IF(.NOT. read_from_restart_file) THEN
        
        CALL init_state()
        
        CALL MPI_BARRIER(comm3d, ierr)

        ! Exchange boundary data with neighbors
        IF(myrank .EQ. 0) PRINT*,''
        IF(myrank .EQ. 0) PRINT*,'Exchanging boundary for initial state.'
     
        CALL start_halo_exchange()
        CALL end_halo_exchange()

        IF(myrank .EQ. 0) PRINT*,'Boundary exchange completed.'
        IF(myrank .EQ. 0) PRINT*,''

        CALL boundary_conditions()
        
        s1_in = s1_def
        s2_in = s2_def
        s3_in = s3_def
        
        dump_count = 0
        
    ELSE
        
        CALL read_restart(t_sim)
        
        s1_in = s1
        s2_in = s2
        s3_in = s3
        
    END IF

    ! initialize turbulence driver
    IF(drive_turbulence) CALL init_turb_driver(s1_in, s2_in, s3_in)
    
    ! compute initial random force 
    IF(drive_turbulence) CALL generate_rand_force()
    
    drive_next = .FALSE.
    
    ! compute time for next turbulence driving
    t_turb = t_sim + turb_frequency
    
    IF(myrank .EQ. 0) PRINT*,'Writing initial state to file...'
    
    CALL MPI_BARRIER(comm3d, ierr)

    ! dump initial state to file
    IF(.NOT. parallel_io) THEN
        CALL writetofile_unformatted(0)
    ELSE
        CALL load_parallel_file_buffer(t_sim)      
        CALL writetofile_unformatted_parallel(0, parallel_io_req, parallel_file_handle, .TRUE., .FALSE., write_status)
    END IF    
    
    IF(myrank .EQ. 0) PRINT*,'Done writing initial state to file...'


    ! initialize time step counter
    tsteps = 1
    
    ! compute times for the next output and restart file dumps
    t_dump = t_sim + dump_frequency
    t_restart = t_sim + restart_frequency 
    dump_count = dump_count + 1

    CALL compute_divB()
   

END SUBROUTINE main_init


! Top-level routine for updating MHD state variables over a time step 
SUBROUTINE mhd_update()

    IF(mhd_solver_type .EQ. 1) THEN
        CALL tvd_solve_3d(my_dt)
    ELSE IF(mhd_solver_type .EQ. 2) THEN
        CALL vanleer_solve_3d(my_dt, drive_next)
    ELSE
        PRINT*,'ERROR! Need to set mhd_solver_type = 1 OR 2'
        STOP
    END IF
 
END SUBROUTINE mhd_update


! Top-level routine for performing a file-dump
SUBROUTINE file_dump(dump_type)

    INTEGER, INTENT(IN) :: dump_type ! 1: state dump, 2: restart dump
    LOGICAL :: write_status

    IF(dump_type .EQ. 1) THEN
    
        IF(print_debug) PRINT*,'Writing MHD state to file...'
            
        IF(.NOT. parallel_io) THEN
            CALL writetofile_unformatted(dump_count)
        ELSE
            CALL load_parallel_file_buffer(t_sim)
            CALL writetofile_unformatted_parallel(dump_count, parallel_io_req, parallel_file_handle, .FALSE., .FALSE., write_status)
        END IF       
       
        IF(write_status) THEN
            WRITE(UNIT=99,FMT='("Dumped MHD state to file. Simulation Time = ",f20.10,", Dump# =",i10)') t_sim, dump_count
        ELSE
            WRITE(UNIT=99,FMT='("Rank # ", i10 ,". MPI parallel file write failed! Unable to dump MHD state to file. Simulation Time = ", &
                  f20.10,", Dump# =",i10)') myrank,t_sim, dump_count
        END IF
       
        IF(print_debug) PRINT*,'Done writing MHD state to file..'
            
        ! increment dump counter    
        dump_count = dump_count + 1

        ! compute the next dump time
        t_dump = MAX(t_dump + dump_frequency, t_sim + 1d-2*my_dt)

    ELSE
    
        CALL write_restart(t_sim)
        
        WRITE(UNIT=99,FMT='("Wrote to restart file. Simulation Time = ",f20.10)') t_sim
        
        ! compute the next dump time
        t_restart = MAX(t_restart + restart_frequency, t_sim + 1d-2*my_dt)
    
    END IF

END SUBROUTINE file_dump


SUBROUTINE compute_new_time_step(dt_new)

    REAL*8, INTENT(INOUT) :: dt_new
    REAL*8 :: cmax

    IF(print_debug) PRINT*,'Computing time step size.'
    
    IF(mhd_solver_type .EQ. 1) THEN
        CALL get_max_speed_tvd(cmax)
    ELSE IF(mhd_solver_type .EQ. 2) THEN
        CALL get_max_speed_vanleer(cmax)
    END IF
    
    dt_new = COUR * my_dx / cmax 


END SUBROUTINE compute_new_time_step



SUBROUTINE share_dt()

   ! initiate non-blocking MPI all_reduce
   CALL start_allreduce(my_dt, reduced_dt, reduce_req, 1)

END SUBROUTINE share_dt


SUBROUTINE get_new_dt() 

   ! wait for non-blocking MPI all_reduce to complete
   CALL end_allreduce(reduce_req) 
  
   ! replace current time step with the global minimum 
   my_dt = MIN(t_end-t_sim, reduced_dt)       
      
END SUBROUTINE get_new_dt


SUBROUTINE compute_divB()

    INTEGER :: i,j,k

    IF(myrank .EQ. 0 .AND. print_debug) PRINT*,'Computing divB..'

    my_divB = 0.d0
    
    DO k = 1, nz
        DO j = 1, ny
            DO i = 1, nx

                my_divB = my_divB + (bface_3d(i,j,k,1) - bface_3d(i-1,j,k,1)) +  (bface_3d(i,j,k,2) - bface_3d(i,j-1,k,2)) + &
                           (bface_3d(i,j,k,3) - bface_3d(i,j,k-1,3))

            END DO        
        END DO        
    END DO        

    my_divB = my_divB / my_dx

    IF(myrank .EQ. 0) PRINT*,'DivB =',my_divB

END SUBROUTINE compute_divB


SUBROUTINE compute_KE()

    INTEGER :: i,j,k
    REAL*8 :: KE
    
    IF(myrank .EQ. 0 .AND. print_debug) PRINT*,'Computing total kinetic energy..'

    KE = 0.d0
    
    DO k = 1, nz
        DO j = 1, ny
            DO i = 1, nx
                KE = KE + (q_3d(i,j,k,2)**2 + q_3d(i,j,k,3)**2 + q_3d(i,j,k,4)**2)/q_3d(i,j,k,1)
            END DO        
        END DO        
    END DO        

    KE = 0.5d0*KE*(my_dx**3)
 
    IF(myrank .EQ. 0) PRINT*,'Total kinetic energy =',KE

END SUBROUTINE compute_KE



SUBROUTINE get_timing_stats()

    REAL*8 :: all_tot, all_tsolve, all_tcom, all_tfile, t_per_step, all_divB 
    INTEGER :: ierr
    CHARACTER(LEN=180) :: filename

    IF(myrank .EQ. 0) PRINT*,'Getting more timing stats...'

    CALL MPI_ALLREDUCE(ttot, all_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
    CALL MPI_ALLREDUCE(tsolve, all_tsolve, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
    CALL MPI_ALLREDUCE(tcom, all_tcom, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
    CALL MPI_ALLREDUCE(tfile, all_tfile, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
    CALL MPI_ALLREDUCE(my_divB, all_divB, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)

    IF(myrank .EQ. 0) THEN
    
        PRINT*,''
        PRINT*,'Average per rank:'
    
        WRITE(*, FMT = '("Total execution time (s) = ",f10.5,", Solve time fraction = ",f8.7,", Total MPI comm. time fraction = ",f8.7, &
              ", Total I/O time fraction = ",f8.7)') all_tot/np, all_tsolve/all_tot, all_tcom/all_tot, all_tfile/all_tot
    
        PRINT*,''
        PRINT*,'Total divB = ',all_divB
        PRINT*,''
        
        filename = TRIM(output_filepath)//TRIM('/Run_Stats/timing_stats.txt')
        
        OPEN(UNIT = 10, FILE = filename)
        
        WRITE(10, FMT= '("Total execution time (s) = ",f10.2)') all_tot/np
        WRITE(10, FMT= '("Solve Time Fraction = ",f15.14)') all_tsolve/all_tot
        WRITE(10, FMT= '("Halo Exchange Time Fraction = ",f15.14)') all_tcom/all_tot
        WRITE(10, FMT= '("I/O Time Fraction = ",f15.14)') all_tfile/all_tot
        WRITE(10, FMT= '("Total divB = ",f15.14)') all_divB
               
        CLOSE(UNIT = 10)
        
    END IF

END SUBROUTINE get_timing_stats



END PROGRAM isothermal_tvd_driver