! 2nd order Van Leer Scheme for 3d passive variable advection (a.k.a. MUSCL-HANCOCK)
! References: Stone 'et al., New Astr., 14, 139-148, 2009  (S09)
!
!-------------------------------------------------------------------------------------------------------------
MODULE passive_solver_mod

USE constants_mod
USE grid_data_mod

IMPLICIT NONE


!*****************************************************************************************************************************!
!                                                     LOCAL VARIABLES                                                         !
!*****************************************************************************************************************************!


! Solver Parameters
REAL*8, PARAMETER  :: ZERO = 1.D-20           ! density floor
REAL*8, PARAMETER :: limiter_beta = 1.D0      ! Sweby limiter: 1 <= beta <= 2

! Solver Local variables 
REAL*8, ALLOCATABLE :: qpass_int(:,:,:),          & ! state-vector 2D work arrays for interface-states
                       w(:,:),                    & ! primitive variables 2D work arrays
                       qpass_work_3d(:,:,:,:),    & ! state-vector 3D work array
                       Fpass(:,:,:,:,:),          & ! time-averaged flux vector
                       Fpass_2d(:,:),             & ! flux buffer for rotations
                       vel_3d(:,:,:,:),           & ! work array for storing velocity
                       vel_tavg_3d(:,:,:,:),      & 
                       vel_2d(:,:)

REAL*8 :: dt, dx

INTEGER :: nmax, nb_pass

!*****************************************************************************************************************************!
!                                                     SUBROUTINES                                                             !
!*****************************************************************************************************************************!

CONTAINS



SUBROUTINE init_passive_solver()
    
    ! make sure we have at least 3 boundary cells
    IF(nb .LT. 3) THEN
        PRINT*,'ERROR! Need nb_pass >= 3 for passive variable solver.'
        STOP
    END IF
    
    nb_pass = 3
    
    CALL allocate_local_variables()


END SUBROUTINE init_passive_solver


SUBROUTINE destroy_passive_solver()

   CALL deallocate_local_variables()
   
END SUBROUTINE destroy_passive_solver


! allocate local variables/work arrays
SUBROUTINE allocate_local_variables()

    INTEGER :: total_grid_memory

    nmax = MAX(nx*nranks_x,ny*nranks_y,nz*nranks_z)

    ! set cell size  
    dx = 1.d0 / DBLE(nmax)
    
    ALLOCATE(qpass_int(1-nb_pass:nmax+nb_pass,1-nb_pass:nmax+nb_pass,2))
    ALLOCATE(w(1-nb_pass:nmax+nb_pass,1-nb_pass:nmax+nb_pass))
    ALLOCATE(qpass_work_3d(1-nb_pass:nx+nb_pass,1-nb_pass:ny+nb_pass,1-nb_pass:nz+nb_pass,npass))
    ALLOCATE(Fpass(1-nb_pass:nx+nb_pass,1-nb_pass:ny+nb_pass,1-nb_pass:nz+nb_pass,npass,3))
    ALLOCATE(Fpass_2d(1-nb_pass:nmax+nb_pass,1-nb_pass:nmax+nb_pass))
    ALLOCATE(vel_3d(1-nb_pass:nx+nb_pass,1-nb_pass:ny+nb_pass,1-nb_pass:nz+nb_pass,3))
    ALLOCATE(vel_tavg_3d(1-nb_pass:nx+nb_pass,1-nb_pass:ny+nb_pass,1-nb_pass:nz+nb_pass,3))
    ALLOCATE(vel_2d(1-nb_pass:nmax+nb_pass,1-nb_pass:nmax+nb_pass))


    total_grid_memory = SIZEOF(qpass_int) + SIZEOF(qpass_work_3d) + SIZEOF(w) + &
                        SIZEOF(Fpass) + SIZEOF(Fpass_2d) + SIZEOF(vel_3d) + &
                        SIZEOF(vel_tavg_3d) + SIZEOF(vel_2d)

    IF(myrank .EQ. 0) THEN
    
    PRINT*,''
    PRINT*,'#############################################################################'
    PRINT*,'Total Memory for passive solver work array storage per MPI rank (Mb) = ',total_grid_memory*1e-6                        
    PRINT*,'#############################################################################'
    PRINT*,''

    END IF
    

END SUBROUTINE allocate_local_variables



SUBROUTINE deallocate_local_variables()

    DEALLOCATE(qpass_int, w, qpass_work_3d)
    DEALLOCATE(Fpass, Fpass_2d, vel_3d, vel_tavg_3d, vel_2d)
 
    
END SUBROUTINE deallocate_local_variables


!*****************************************************************************************************************************!

! Top-level Isothermal TVD MHD solve routine (using CTU + CT scheme)
SUBROUTINE passive_solve_3d(dt_in)

    REAL*8, INTENT(INOUT) :: dt_in
    REAL*8 :: dtdx
    INTEGER :: i, j, offset  
    
    dt = dt_in
    
    ! First save initial state in work arrays
    qpass_work_3d = qpass_3d

    offset = 0
    
    !*****************************************************************************************************************
    ! Step 1: Compute first order predictor cell-interafce fluxes along x, y and z directions using 1D Riemann solver
    !*****************************************************************************************************************
    
    IF(print_debug) PRINT*,'Computing predictor TVD x-fluxes.'
    
    CALL compute_xfluxes_3d(1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset, .FALSE.) 

    IF(print_debug) PRINT*,'Computing predictor TVD y-fluxes.'
    
    CALL compute_yfluxes_3d(1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset, .FALSE.) 

    IF(print_debug) PRINT*,'Computing predictor TVD z-fluxes.'
    
    CALL compute_zfluxes_3d(1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset, .FALSE.) 
    
   

    !******************************************************************
    ! Step 3: Apply 1/2 dt update to cell-center passive variables
    !******************************************************************

    dtdx = 0.5d0 * dt / dx
   
    IF(print_debug) PRINT*,'Computing hydro state half update.'
    
    CALL update_passive_state_3d(dtdx, vel_3d, 1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset)

    offset = 1
    
    !************************************************************************************************
    ! Step 5 & 6: Apply piecewise linear reconstruction to compute interface-states at half time-step
    !             then use these to compute the corrector interface fluxes
    !********************************************************************************************
       
    ! compute corrector x-fluxes    
    IF(print_debug) PRINT*,'Computing corrector fluxes in y direction.'
    
    CALL compute_xfluxes_3d(1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset, .TRUE.)  
    
     
    ! compute corrector y-fluxes     
    IF(print_debug) PRINT*,'Computing corrector fluxes in y direction.'

    CALL compute_yfluxes_3d(1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset, .TRUE.)      

  
    ! compute corrector z-fluxes
    IF(print_debug) PRINT*,'Computing corrector fluxes in z direction.'

    CALL compute_zfluxes_3d(1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset, .TRUE.) 

    offset = 2

    !***********************************************************************************
    ! Step 9: Apply final full-dt update to cell-center passive variables
    !*********************************************************************************** 

    ! Restore intial state 
    qpass_3d = qpass_work_3d
        
    dtdx = dt / dx

    IF(print_debug) PRINT*,'Computing final unsplit update of passive variables.'

    CALL update_passive_state_3d(dtdx, vel_tavg_3d, 1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset)

    IF(print_debug) PRINT*,'Update for time step completed.'
 
    
END SUBROUTINE passive_solve_3d


SUBROUTINE copy_velocity(predictor_stage)

    LOGICAL, INTENT(IN) :: predictor_stage
    INTEGER :: i, j, k
    REAL*8 :: irho

    IF(.NOT. predictor_stage) THEN
        
        DO k = 1-nb_pass, nz+nb_pass 
            DO j = 1-nb_pass, ny+nb_pass 
                DO i = 1-nb_pass, nx+nb_pass 
         
                    irho = 1.d0 / q_3d(i,j,k,1)
                    vel_3d(i,j,k,1) = q_3d(i,j,k,2) * irho 
                    vel_3d(i,j,k,2) = q_3d(i,j,k,3) * irho 
                    vel_3d(i,j,k,3) = q_3d(i,j,k,4) * irho 
        
                END DO 
            END DO 
        END DO 
        
    ELSE
    
      ! compute time-averaged velocity for predictor stage
      DO k = 1-nb_pass, nz+nb_pass 
            DO j = 1-nb_pass, ny+nb_pass 
                DO i = 1-nb_pass, nx+nb_pass 
         
                    irho = 1.d0 / q_3d(i,j,k,1)
                    vel_tavg_3d(i,j,k,1) = 0.5d0 * (vel_3d(i,j,k,1) + q_3d(i,j,k,2) * irho) 
                    vel_tavg_3d(i,j,k,2) = 0.5d0 * (vel_3d(i,j,k,2) + q_3d(i,j,k,3) * irho)  
                    vel_tavg_3d(i,j,k,3) = 0.5d0 * (vel_3d(i,j,k,3) + q_3d(i,j,k,4) * irho)
        
                END DO 
            END DO 
        END DO 
    
    END IF
    

END SUBROUTINE copy_velocity



SUBROUTINE compute_xfluxes_3d(ilow, ihi, jlow, jhi, klow, khi, predictor_stage) 

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi, klow, khi
    LOGICAL, INTENT(IN) :: predictor_stage
    INTEGER :: i, j, k, ipass, offset 
    
    ! loop over passive variables
    DO ipass = 1, npass
 
        ! Loop over z planes    
        DO k = klow-nb_pass, khi+nb_pass
     
            IF(print_debug) PRINT*,' Z plane index,  k = ', k
            
            ! Compute interface states. (For predictor stage, the interface states needed for first order fluxes
            ! are just cell-center values) 
            IF(.NOT. predictor_stage) THEN
                
                DO j = jlow-nb_pass, jhi+nb_pass 
                    DO i = ilow-nb_pass, ihi+nb_pass-1        

                        ! left state
                        qpass_int(i,j,1) = qpass_3d(i,j,k,ipass)

                        ! right state
                        qpass_int(i,j,2) = qpass_3d(i+1,j,k,ipass)
                        
                    END DO
                END DO
                
                ! copy velocity at beginning of time step into work array
                DO j = jlow-nb_pass, jhi+nb_pass 
                    DO i = ilow-nb_pass, ihi+nb_pass        

                        vel_2d(i,j) = vel_3d(i,j,k,1) 
                        
                    END DO
                END DO
                                
                offset = 0 

            ELSE
            
                ! compute cell-center value
                DO j = jlow-nb_pass, jhi+nb_pass 
                    DO i = ilow-nb_pass, ihi+nb_pass                     
                              
                        w(i,j) = qpass_3d(i,j,k,ipass)
                
                    END DO
                END DO
                                
                ! Now interpolate to get cell-iterface primitive variables and reconstruct cell-interface states from that 
                ! (for corrector stage, generating second order fluxes requires that we linearly interpolate from cell-center to the interface values)
                CALL compute_interface_states(ilow, ihi, jlow, jhi)
                

                ! copy velocity at half time step into work array
                DO j = jlow-nb_pass, jhi+nb_pass 
                    DO i = ilow-nb_pass, ihi+nb_pass        

                        vel_2d(i,j) = vel_tavg_3d(i,j,k,1) 
                        
                    END DO
                END DO
                
                
                offset = 1
            
            END IF


            CALL compute_fluxes_1d(Fpass_2d, vel_2d, ilow+offset, ihi-offset, jlow, jhi)
            
            
            DO j = jlow-nb_pass, jhi+nb_pass 
                DO i = ilow+offset-nb_pass, ihi-offset+nb_pass-1   
                    Fpass(i,j,k,ipass,1) = Fpass_2d(i,j)
                END DO
            END DO
                
        
        END DO
        
    END DO

END SUBROUTINE compute_xfluxes_3d 


SUBROUTINE compute_yfluxes_3d(ilow, ihi, jlow, jhi, klow, khi, predictor_stage) 

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi, klow, khi
    LOGICAL, INTENT(IN) :: predictor_stage
    INTEGER :: i, j, k, ipass, offset

    ! loop over passive variables
    DO ipass = 1, npass
 
        ! loop over z planes
        DO k = klow-nb_pass, khi+nb_pass
        
            IF(print_debug) PRINT*,' Z plane index,  k = ', k
       

            ! Compute interface states. (For predictor stage, the interface states needed for first order fluxes
            ! are just cell-center values)      
            IF(.NOT. predictor_stage) THEN
                
                ! rotate the stave vector array so that fastest index runs along y-direction                
                DO i = ilow-nb_pass, ihi+nb_pass 
                    DO j = jlow-nb_pass, jhi+nb_pass-1
                
                        ! left state
                        qpass_int(j,i,1) = qpass_3d(i,j,k,ipass)
                       
                        ! right state
                        qpass_int(j,i,2) = qpass_3d(i,j+1,k,ipass)
                       
                    END DO
                END DO
         
                ! copy velocity at beginning of time step into work array
                DO i = ilow-nb_pass, ihi+nb_pass 
                    DO j = jlow-nb_pass, jhi+nb_pass       

                        vel_2d(j,i) = vel_3d(i,j,k,2) 
                            
                    END DO
                END DO
                                    
                offset = 0
         
            ELSE
            
                ! compute cell-center values
                ! also rotate so that fastest index runs along y-direction
                DO i = ilow-nb_pass, ihi+nb_pass 
                    DO j = jlow-nb_pass, jhi+nb_pass                                       
                        w(j,i) = qpass_3d(i,j,k,ipass)                 
                    END DO
                END DO
                
                ! Now interpolate to get cell-iterface primitive variables and reconstruct cell-interface states from that 
                ! (for corrector stage, generating second order fluxes requires that we linearly interpolate from cell-center to the interface values)
                CALL compute_interface_states(jlow, jhi, ilow, ihi)
                    
         
                ! copy velocity at half time step into work array
                DO i = ilow-nb_pass, ihi+nb_pass 
                    DO j = jlow-nb_pass, jhi+nb_pass       

                        vel_2d(j,i) = vel_tavg_3d(i,j,k,2) 
                            
                    END DO
                END DO
                
                offset = 1      
            
            END IF
         
            
            CALL compute_fluxes_1d(Fpass_2d, vel_2d, jlow+offset, jhi-offset, ilow, ihi)
           
           
            ! Rotate the y-fluxes in flux buffer and store them in main flux array 
            DO j = jlow+offset-nb_pass, jhi-offset+nb_pass-1
                DO i = ilow-nb_pass, ihi+nb_pass
                    Fpass(i,j,k,ipass,2) = Fpass_2d(j,i)             
                END DO
            END DO
            
        END DO

    END DO
    

END SUBROUTINE compute_yfluxes_3d 


SUBROUTINE compute_zfluxes_3d(ilow, ihi, jlow, jhi, klow, khi, predictor_stage) 

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi, klow, khi
    LOGICAL, INTENT(IN) :: predictor_stage
    INTEGER :: i, j, k, ipass, offset

    ! loop over passive variables
    DO ipass = 1, npass

        ! loop over y planes
        DO j = jlow-nb_pass, jhi+nb_pass  
        
            IF(print_debug) PRINT*,' Y plane index:  j = ', j
        
            ! Compute interface states. (For predictor stage, the interface states needed for first order fluxes
            ! are just cell-center values)      
            IF(.NOT. predictor_stage) THEN
                
                
                ! first rotate the stave vector array so that fastest index runs along z-direction 
                DO i = ilow-nb_pass, ihi+nb_pass 
                    DO k = klow-nb_pass, khi+nb_pass-1
                
                        ! left state
                        qpass_int(k,i,1) = qpass_3d(i,j,k,ipass)
                
                        ! right state
                        qpass_int(k,i,2) = qpass_3d(i,j,k+1,ipass)
                
                    END DO
                END DO
                
                ! copy velocity at beginning of time step into work array
                DO i = ilow-nb_pass, ihi+nb_pass 
                    DO k = klow-nb_pass, khi+nb_pass
                        vel_2d(k,i) = vel_3d(i,j,k,3)   
                    END DO
                END DO
                
                offset = 0
                
            ELSE
            
                ! compute cell-center values
                DO i = ilow-nb_pass, ihi+nb_pass 
                    DO k = klow-nb_pass, khi+nb_pass                                                     
                        w(k,i) = qpass_3d(i,j,k,ipass)
                    END DO
                END DO
                
                ! Now interpolate to get cell-iterface primitive variables and reconstruct cell-interface states from that 
                ! (for corrector stage, generating second order fluxes requires that we linearly interpolate from cell-center to the interface values)
                CALL compute_interface_states(klow, khi, ilow, ihi)
                        
                ! copy velocity at half time step into work array
                DO i = ilow-nb_pass, ihi+nb_pass 
                    DO k = klow-nb_pass, khi+nb_pass
                        vel_2d(k,i) = vel_tavg_3d(i,j,k,3)   
                    END DO
                END DO
                
                offset = 1      
            
            END IF
           
            CALL compute_fluxes_1d(Fpass_2d, vel_2d, klow+offset, khi-offset, ilow, ihi)
            
            ! Rotate the z-fluxes in flux buffer and store them in main flux array 
            DO k = klow+offset-nb_pass, khi-offset+nb_pass-1
                DO i = ilow-nb_pass, ihi+nb_pass        
                    Fpass(i,j,k,ipass,3) = Fpass_2d(k,i)               
                END DO
            END DO

        END DO

    END DO


END SUBROUTINE compute_zfluxes_3d 



SUBROUTINE compute_interface_states(ilow, ihi, jlow, jhi)

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi
    INTEGER :: i, j, k
    REAL*8 :: delw_L, delw_R, delw_C, delw(1-nb_pass:nmax+nb_pass,1-nb_pass:nmax+nb_pass)
    REAL*8 :: signL, signR, signC
    
    ! Compute (cell-center) left, right and centered differences of the passive variable and the monotonized slope 
    DO j = jlow-nb_pass, jhi+nb_pass
        DO i = ilow-nb_pass+1, ihi+nb_pass-1

            ! compute left difference
            delw_L = w(i,j) - w(i-1,j) 
            
            ! compute right difference
            delw_R = w(i+1,j) - w(i,j) 
    
            ! compute centered difference
            !delw_C = 0.5d0 *( delw_L + delw_R )
    
            ! compute monotonized slope  (min-mod limiter)
            delw(i,j) =  0.5d0 * (SIGN(1.d0, delw_L) + SIGN(1.d0, delw_R)) * MIN(ABS(delw_L), ABS(delw_R))                         
            
        END DO
    END DO

    ! compute left and right interface states
    DO j = jlow-nb_pass, jhi+nb_pass
        DO i = ilow-nb_pass+1, ihi+nb_pass-2

            ! left state
            qpass_int(i,j,1) =  w(i,j) + 0.5d0 * delw(i,j)  

        END DO
    END DO

    DO j = jlow-nb_pass, jhi+nb_pass
        DO i = ilow-nb_pass+1, ihi+nb_pass-2

            ! right state
            qpass_int(i,j,2) =  w(i+1,j) - 0.5d0 * delw(i+1,j)  

        END DO
    END DO
    
    
    
END SUBROUTINE compute_interface_states



! (unsplit) finite volume update of the passive variable 
SUBROUTINE update_passive_state_3d(dtdx, vel, ilow, ihi, jlow, jhi, klow, khi)

    REAL*8, INTENT(IN) :: dtdx
    REAL*8, INTENT(IN) :: vel(1-nb_pass:nx+nb_pass,1-nb_pass:ny+nb_pass,1-nb_pass:nz+nb_pass,3)
    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi, klow, khi
    INTEGER :: i, j, k, ipass
  

    ! loop over passive variables
    DO ipass = 1, npass

        ! update the conserved variables (density and momentum only)
         DO k = klow-nb_pass+1, khi+nb_pass-1    
             DO j = jlow-nb_pass+1, jhi+nb_pass-1     
                 DO i = ilow-nb_pass+1, ihi+nb_pass-1
            
                    qpass_3d(i,j,k,ipass) = qpass_3d(i,j,k,ipass) - dtdx *  ( &
                                            vel(i,j,k,1) * (Fpass(i,j,k,ipass,1) - Fpass(i-1,j,k,ipass,1)) + &
                                            vel(i,j,k,2) * (Fpass(i,j,k,ipass,2) - Fpass(i,j-1,k,ipass,2)) + &
                                            vel(i,j,k,3) * (Fpass(i,j,k,ipass,3) - Fpass(i,j,k-1,ipass,3)) ) 
                                    

                END DO    
            END DO
        END DO

    END DO

END SUBROUTINE update_passive_state_3d


! This subroutine computes 1D upwind flux
SUBROUTINE compute_fluxes_1d(Flux, vx, ilow, ihi, jlow, jhi)

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi
    REAL*8, INTENT(INOUT) :: Flux(1-nb_pass:nmax+nb_pass,1-nb_pass:nmax+nb_pass)
    REAL*8, INTENT(IN) :: vx(1-nb_pass:nmax+nb_pass,1-nb_pass:nmax+nb_pass)
    INTEGER :: i, j, k
    REAL*8 :: uL, uR, vavg, theta
  
  
    DO j = jlow-nb_pass, jhi+nb_pass
        DO i = ilow-nb_pass, ihi+nb_pass-1
        
     
            ! left and right states
            uL = qpass_int(i,j,1)
            uR = qpass_int(i,j,2)
            
            ! fluid velocity at interface
            vavg = 0.5d0 * (vx(i,j) + vx(i+1,j))
      
            IF(vavg .GE. 0.d0) THEN
                theta = 1.d0                    
            ELSE
                theta = -1.d0
            END IF
      
            ! compute upwind flux
            Flux(i,j) = 0.5d0 * ( (1.d0 + theta) * uL + (1.d0 - theta) * uR )                    

        END DO    
    END DO   
    
    
END SUBROUTINE compute_fluxes_1d



END MODULE passive_solver_mod