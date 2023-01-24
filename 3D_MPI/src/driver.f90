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
INTEGER :: KE_sum_req, vrms_sum_req, dt_req, mx_req, my_req, mz_req, parallel_io_req, parallel_file_handle
REAL*8 ::  COUR, t_old, t_dump, t_restart, t_turb   
LOGICAL :: write_status, drive_next
REAL*8 :: my_KE, my_vrms, my_divB, my_momx, my_momy, my_momz, sum_momx, sum_momy, sum_momz, sum_KE, sum_vrms 

!########################################################################################################################
   
! do initializations for grid, solver, mpi and i/o to get everything started   
CALL main_init()
    
CALL MPI_BARRIER(comm3d, ierr)

IF(myrank .EQ. 0) PRINT*,'Starting simulation loop.'

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
    
    ! set courant number to default value
    COUR = COUR_def   
   
    ! evolve MHD state vector and passive variables
    tsolve_1 = MPI_Wtime()
    IF(myrank .EQ. 0) PRINT*,'MHD solve in progress.'

    IF(npass .GT. 0) CALL copy_velocity(.FALSE.)
    
    IF(drive_turbulence .AND. (t_turb .GT. t_old) .AND. (t_turb .LE. t_old+my_dt)) THEN
        drive_next = .TRUE.
        IF(myrank .EQ. 0) THEN
            PRINT*,'##############################################'        
            PRINT*, 'Applying turbulence driving force.'
            PRINT*,'##############################################'
        END IF

        ! compute the next turbulence driving time
        t_turb = MAX(t_turb + dt_turb, t_sim + 1d-2*my_dt)
           
        my_dt = 0.5d0 * my_dt ! maybe unnecessary, but just to be safe...
           
    END IF   
        
    CALL mhd_update()
    
    IF(npass .GT. 0) THEN  !##### Warning: With trubulence forcing enabled, need to figure out whether to exclude the force from the time averaged velocity used for advecting passive vars
        CALL copy_velocity(.TRUE.)
        IF(myrank .EQ. 0) PRINT*,'Updating passive variables'
        CALL passive_solve_3d(my_dt)
        IF(myrank .EQ. 0) PRINT*,'Done updating passive variables'
    END IF        

    IF(myrank .EQ. 0) PRINT*,'MHD solve completed.'
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
    tsolve_3 = MPI_Wtime()
    IF(drive_turbulence .AND. drive_next) THEN
        IF(myrank .EQ. 0) PRINT*,'Computing random force for next drive cycle...'
        CALL generate_rand_force()
        IF(myrank .EQ. 0) PRINT*,'Done computing random force.'
        drive_next = .FALSE.    
        
        ! lower the courant number for next time step to avoid generation of spurious nunmerical oscillations due to forcing
        COUR = 0.5d0*COUR_def  
        
    END IF    
    tsolve_4 = MPI_Wtime()
    tsolve = tsolve + tsolve_4 - tsolve_3
    
    ! save state to file 
    tfile_1 = MPI_Wtime()
    write_status = .FALSE.
    IF((t_dump .GT. t_old) .AND. (t_dump .LE. t_sim)) THEN
        CALL file_dump(1)
    END IF
    tfile_2 = MPI_Wtime()
    
    tfile_3 = MPI_Wtime()
    IF((t_restart .GT. t_old) .AND. (t_restart .LE. t_sim)) THEN
        CALL file_dump(2)
    END IF
    tfile_4 = MPI_Wtime() 
    
    ! compute new local time step
    tsolve_3 = MPI_Wtime()
    CALL compute_new_time_step(my_dt)    
    
    ! some diagnostics    
    CALL compute_KE()
    !CALL compute_divB() 
    tsolve_4 = MPI_Wtime()
    tsolve = tsolve + tsolve_4 - tsolve_3

  
    
    !*****************************************************************************************************
    !*****************************************************************************************************
    tcom_3 = MPI_Wtime()
    CALL end_halo_exchange()
    IF(myrank .EQ. 0) PRINT*,'Halo exchange completed.'
    tcom_4 = MPI_Wtime()

        
    ! post non-blocking mpi all reduce to communicate my time step value with all other ranks    
    IF(myrank .EQ. 0) PRINT*,'Initiating mpi_all_reduce.'

    CALL share_dt() 
    
    !*****************************************************************************************************
    ! do some other work here instead of just waiting for allreduce to complete...
    !*****************************************************************************************************
    !*****************************************************************************************************

    ! apply boundary conditions (make sure to do this after halo exchange with neighbor's is done)
    CALL boundary_conditions()
        
    !*****************************************************************************************************
    !*****************************************************************************************************
    ! wait for MPI reduction to complete and get the reduced time step value
    IF(myrank .EQ. 0) PRINT*,'Waiting for all_reduce to complete...'
    tcom_5 = MPI_Wtime()
    CALL finish_mpi_reductions()
    tcom_6 = MPI_Wtime()    

    IF(myrank .EQ. 0) PRINT*,'Obtained new time step from MPI all reduce.'   

    tfile = tfile + tfile_2-tfile_1 + tfile_4-tfile_3    
    trma = trma + tcom_2-tcom_1 + tcom_4-tcom_3
    tred = tred + tcom_4-tcom_3 + tcom_6-tcom_5

    IF(myrank .EQ. 0) THEN
        PRINT*,'Step Total Time (sec) = ',tcom_6-t2
        PRINT*,'Step Solve Time (sec) = ',tsolve_2-tsolve_1 + tsolve_4-tsolve_3
        PRINT*,'Step Halo Exchange Time (sec) = ',tcom_4-tcom_1
        PRINT*,'Step ALLREDUCE Wait Time (sec) = ',tcom_6-tcom_4
        PRINT*,'Step MHD State File I/O Time (sec) = ',tfile_2-tfile_1
        PRINT*,'Step Restart File I/O Time (sec) = ',tfile_4-tfile_3
        PRINT*,''
    END IF

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


IF(myrank .EQ. 0) PRINT*,'Done!'


CONTAINS


! Top-level routine for initialization of every major component 
SUBROUTINE main_init()

    INTEGER :: s1_in, s2_in, s3_in, ierr


    ! start MPI and set up domain decomposition
    CALL initialize_MPI()

    ! allocate memory for grid data
    CALL create_grid_arrays()

    ! set initial time step size
    my_dt = 1.d-15 
    t_sim = 0.d0
    COUR = COUR_def
    
    ! set cell size
    my_dx = 1.d0 / DBLE(MAX(nx*nranks_x,ny*nranks_y,nz*nranks_z))


    IF(myrank .EQ. 0) PRINT*,'Setting up grid and MPI domain decomposition.'
 
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

        IF(myrank .EQ. 0) PRINT*,'Applying boundary conditions'
        CALL boundary_conditions()
        
        s1_in = s1_def
        s2_in = s2_def
        s3_in = s3_def
        
        dump_count = 0
        
    ELSE
        
        CALL read_restart(t_sim)
        
        CALL MPI_BARRIER(comm3d, ierr)

        ! Exchange boundary data with neighbors
        IF(myrank .EQ. 0) PRINT*,''
        IF(myrank .EQ. 0) PRINT*,'Exchanging boundary for initial state.'

     
        CALL start_halo_exchange()
        CALL end_halo_exchange()

        IF(myrank .EQ. 0) PRINT*,'Boundary exchange completed.'
        IF(myrank .EQ. 0) PRINT*,''

        IF(myrank .EQ. 0) PRINT*,'Applying boundary conditions'
        CALL boundary_conditions()
        
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
    t_turb = t_sim + dt_turb
    
    IF(myrank .EQ. 0) PRINT*,'Writing initial state to file...'
    
    CALL MPI_BARRIER(comm3d, ierr)


    ! dump initial state to file
    IF(.NOT. parallel_io) THEN
        CALL writetofile_unformatted(dump_count)
    ELSE
        CALL load_parallel_file_buffer(t_sim)      
        CALL writetofile_unformatted_parallel(dump_count, parallel_io_req, parallel_file_handle, .TRUE., .FALSE., write_status)
    END IF    
    
    
    IF(myrank .EQ. 0) PRINT*,'Done writing initial state to file...'

    ! initialize time step counter
    tsteps = 1
    
    ! compute times for the next output and restart file dumps
    t_dump = t_sim + dump_frequency
    t_restart = t_sim + restart_frequency 
    IF(.NOT. read_from_restart_file) dump_count = dump_count + 1

    CALL compute_divB()
   
    IF(myrank .EQ. 0) CALL output_runpars()

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

    IF(dump_type .EQ. 1) THEN
    
        IF(myrank .EQ. 0) PRINT*,'Dumping MHD state to output file...'
            
        IF(.NOT. parallel_io) THEN
            CALL writetofile_unformatted(dump_count)
            write_status = .TRUE.            
        ELSE
            CALL load_parallel_file_buffer(t_sim)
            CALL writetofile_unformatted_parallel(dump_count, parallel_io_req, parallel_file_handle, .FALSE., .FALSE., write_status)
        END IF       
       
       
        IF(myrank .EQ. 0) PRINT*,'Done dumping MHD state to file..'
            
        ! increment dump counter    
        dump_count = dump_count + 1

        ! compute the next dump time
        t_dump = MAX(t_dump + dump_frequency, t_sim + 1d-2*my_dt)

    ELSE
    
        CALL write_restart(t_sim)
                
        ! compute the next dump time
        t_restart = MAX(t_restart + restart_frequency, t_sim + 1d-2*my_dt)
    
        IF(myrank .EQ. 0) PRINT*,'Wrote to restart file.'
    
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
   CALL start_allreduce(my_dt, reduced_dt, dt_req, 1)
   CALL start_allreduce(my_vrms, sum_vrms, vrms_sum_req, 2)
   CALL start_allreduce(my_KE, sum_KE, KE_sum_req, 2)
   CALL start_allreduce(my_momx, sum_momx, mx_req, 2)
   CALL start_allreduce(my_momy, sum_momy, my_req, 2)
   CALL start_allreduce(my_momz, sum_momz, mz_req, 2)
   
END SUBROUTINE share_dt


SUBROUTINE finish_mpi_reductions() 

    ! wait for non-blocking MPI all_reduce to complete
    CALL end_allreduce(dt_req) 
    CALL end_allreduce(vrms_sum_req) 
    CALL end_allreduce(KE_sum_req) 
    CALL end_allreduce(mx_req) 
    CALL end_allreduce(my_req) 
    CALL end_allreduce(mz_req) 
  
    ! replace current time step with the global minimum 
    my_dt = MIN(t_end-t_sim, reduced_dt)       
      
    sum_vrms = SQRT(sum_vrms)
    
    IF(myrank .EQ. 0) THEN
        PRINT*,''
        PRINT*,'<v_rms> =',sum_vrms
        PRINT*,'<kinetic energy> =',sum_KE
        PRINT*,'<momentum x> = ', sum_momx
        PRINT*,'<momentum y> = ', sum_momy
        PRINT*,'<momentum z> = ', sum_momz
        PRINT*,''
    END IF
      
END SUBROUTINE finish_mpi_reductions


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

    !my_divB = my_divB / my_dx

    !IF(myrank .EQ. 0) PRINT*,'DivB =',my_divB

END SUBROUTINE compute_divB


SUBROUTINE compute_KE()

    INTEGER :: i,j,k
    REAL*8 :: vsqr, rho
    
    IF(myrank .EQ. 0 .AND. print_debug) PRINT*,'Computing total kinetic energy..'

    my_KE = 0.d0
    my_vrms = 0.d0
    my_momx = 0.d0
    my_momy = 0.d0
    my_momz = 0.d0
    
    DO k = 1, nz
        DO j = 1, ny
            DO i = 1, nx
                rho = MIN(MAX(q_3d(i,j,k,1), 1.d-15),1.d10)
                vsqr = MIN((q_3d(i,j,k,2)**2 + q_3d(i,j,k,3)**2 + q_3d(i,j,k,4)**2)/(rho**2), 1.d5)
                my_vrms = my_vrms + vsqr
                my_KE = my_KE + rho*vsqr
                my_momx = my_momx + q_3d(i,j,k,2)
                my_momy = my_momy + q_3d(i,j,k,3)
                my_momz = my_momz + q_3d(i,j,k,4)
                
            END DO        
        END DO        
    END DO        

    my_vrms = my_vrms * (my_dx**3)
    my_KE = 0.5d0 * my_KE * (my_dx**3)
    my_momx = my_momx * (my_dx**3)
    my_momy = my_momy * (my_dx**3)
    my_momz = my_momz * (my_dx**3)
 
   
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


SUBROUTINE output_runpars()

    CHARACTER(LEN=180) :: filename
    
    filename = TRIM(output_filepath)//TRIM('/run_pars.txt')
        
    OPEN(UNIT=10, FILE = filename)
        
    WRITE(10, FMT= '("Cube Size = ",i6,"x",i6,"x",i6)') nx*nranks_x,ny*nranks_y,nz*nranks_z
    WRITE(10, FMT= '("t_end = ",f8.2)') t_end
    WRITE(10, FMT= '("Isothermal sound speed = ",f8.2)') sound_speed
    WRITE(10, FMT= '("Edot = ",e9.2)') Edot
    WRITE(10, FMT= '("dt_turb = ",f8.4)') dt_turb
    WRITE(10, FMT= '("k_max = ",i2)') k_max
    WRITE(10, FMT= '("Bx_0 = ",e9.2)') Bx_0
    WRITE(10, FMT= '("By_0 = ",e9.2)') By_0
    WRITE(10, FMT= '("Bz_0 = ",e9.2)') Bz_0
    WRITE(10, FMT= '("Output dump frequency = ",f8.4)') dump_frequency
    WRITE(10, FMT= '("Restart frequency = ",f8.4)') restart_frequency
               
    CLOSE(UNIT=10)


END SUBROUTINE output_runpars


END PROGRAM isothermal_tvd_driver