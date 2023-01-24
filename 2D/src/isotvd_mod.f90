! 1D Isothermal MHD 2nd order TVD Solver
! References: Kim 'et al, 1998  (K98)
!             Stone et al, 2018 (S18)
!
!-----------------------------------------
! Outline of 2D CTU+CT Algorithm from S18:
!-----------------------------------------
! 
! The goal is to do a directionally unsplit update of the Isothermal MHD state variables and
! maintain div(B) = 0 to machine precision.
!
! ****************
! * CT Algorithm * 
! ****************
! Div(B) = 0 is maintained by evolving area-averaged transverse magnetic fields 
! defined at each cell-faces via circulation of longitudinal electric fields fields
! at the bounding cell edges, or cell-corners in 2D. (This is just the area integral
! version of Faraday's Law.). Details on how to compute the cell edge/corner electric 
! fields are given in (S18).
!
! 
! *****************
! * CTU Algorithm * 
! *****************
!
! The 2D update at each time step is performed simulatenously in both 
! x,y directions using a (sort of) predictor corrector approach.
!
! 1) Predictor Stage: We first perform directionally split 1/2 time step updates on the 
!                     MHD state variables (except B field). The update in a given 
!                     direction is done using the flux in the transverse direction.
!                     The B field undergoes a half time step CT update.
!                     (Note: The dimensionally split 1/2 time updates require extra "multidimensional
!                     source terms" to be included. These terms are proportional to the derivative of B-longitudinal
!                     along the given direction).
!
! 2) Corrector Stage: The MHD state variables and cell-face magnetic fields obtained from 
!                     the half-update in the predictor stage are then used to compute "corrector"  
!                     TVD fluxes and cell-corner (cell-edges in 3D) electric electric fields.
!                     These "corrector" fluxes and electric fields are then used to fully update 
!                     the MHD state variables without using any dimensional splitting, i.e. the 
!                     x and y TVD flux gradients are simultaneously applied (and the CT magnetic 
!                     field update is already dimensionally unsplit). 
!
!                     The 1D TVD fluxes are computed using the 2nd order Roe-type Riemann solver in K98.
!
!-------------------------------------------------------------------------------------------------------------


MODULE Iso_TVD_Solver_mod

USE constants_mod
USE grid_data_mod

IMPLICIT NONE

! Solver Parameters
INTEGER, PARAMETER :: nwaves = 6
REAL*8, PARAMETER  :: ZERO = 1.D-20 
REAL*8, PARAMETER :: limiter_beta = 1.D0  ! Sweby limiter: 1 <= beta <= 2
REAL*8, PARAMETER :: COUR = 0.7D0 ! Courant number

! Solver Local variables 
REAL*8, ALLOCATABLE :: q(:,:,:),       &      ! state-vector work array
                       bface(:,:,:),   &      ! cell-face magnetic field work array 
                       Fs(:,:,:,:),    &      ! time-averaged flux vector
                       Fs_buff(:,:,:,:),  &   ! flux buffer for rotations
                       qavg(:,:,:),  &        ! Roe-averaged state    
                       Rk(:,:,:,:),  &        ! Roe-averaged right eigenvectors
                       ck(:,:,:),    &        ! characteristic variables (i.e. projection of state vector onto the left eigenvectors)   
                       eigenvalues(:,:,:),  & ! Roe-averaged eigenvalues
                       wave_speeds(:,:,:),  & ! Roe-averaged wave speeds
                       emf_2d(:,:,:),       & ! cell-center reference and cell-face electric field z-component
                       emf_corner(:,:)        ! cell-corner CT electric field z-component
                       
REAL*8, PARAMETER :: eps(6) = (/ 0.3, 0.0, 0.3, 0.3, 0.0 ,0.3/) ! wave dissipation constants

REAL*8 :: dt, dx

INTEGER :: nmax

CONTAINS


SUBROUTINE init_solver()
    
    !OPEN(UNIT = 12, FILE='Output/flux.txt')
    CALL allocate_local_variables()


END SUBROUTINE init_solver


SUBROUTINE destroy_solver()

   CALL deallocate_local_variables()
   !CLOSE(UNIT=12)
   
END SUBROUTINE destroy_solver

! allocate local variables/work arrays
SUBROUTINE allocate_local_variables()

    nmax = MAX(nx,ny)

    ALLOCATE(q(1-nb:nmax+nb,1-nb:nmax+nb,7))
    ALLOCATE(bface(1-nb:nx+nb, 1-nb:ny+nb,2))
    ALLOCATE(Fs(1-nb:nmax+nb,1-nb:nmax+nb,7,2))
    ALLOCATE(Fs_buff(1-nb:nmax+nb,1-nb:nmax+nb,7,2))
    ALLOCATE(qavg(1-nb:nmax+nb,1-nb:nmax+nb,7))
    ALLOCATE(Rk(1-nb:nmax+nb,1-nb:nmax+nb,nwaves, nwaves))
    ALLOCATE(ck(1-nb:nmax+nb,1-nb:nmax+nb,nwaves))
    ALLOCATE(eigenvalues(1-nb:nmax+nb,1-nb:nmax+nb,nwaves))
    ALLOCATE(wave_speeds(1-nb:nmax+nb,1-nb:nmax+nb,3))

    ALLOCATE(emf_2d(1-nb:nx+nb,1-nb:ny+nb,3))
    ALLOCATE(emf_corner(1-nb:nx+nb,1-nb:ny+nb))


END SUBROUTINE allocate_local_variables


SUBROUTINE deallocate_local_variables()

    DEALLOCATE(q, bface, Fs, Fs_buff, qavg, Rk, ck, eigenvalues, wave_speeds, emf_2d, emf_corner)
    
END SUBROUTINE deallocate_local_variables

! Top-level Isothermal TVD MHD solve routine (using CTU + CT scheme)
SUBROUTINE solve_2d()

    REAL*8 :: cmax
    INTEGER :: i, j
    
    
    dx = 1.d0 / nx
    dt = 0.d0


    ! clear electric field arrays
    emf_2d = 0.d0
    emf_corner = 0.d0

    ! copy initial cell face magnetic field into work array
    bface = bface_2d


    !************************
    ! Compute time-step size
    !************************
    PRINT*,'Computing time step size.'
    CALL get_max_speed(cmax)
    
    dt = COUR * dx / cmax 
    PRINT*,'dt =',dt
    PRINT*,'dtdx = ', dt/dx



    !*******************************************************************************
    ! Step 1: Compute 1/2 dt TVD fluxes along x and y directions using 1D TVD solver
    !*******************************************************************************
    PRINT*,'Computing predictor TVD x-fluxes.'
    
    CALL compute_xfluxes_2d(.FALSE.)  ! send false flag since this is the predictor stage

    PRINT*,'Computing predictor TVD y-fluxes.'
    
    CALL compute_yfluxes_2d(.FALSE.)  ! send false flag since this is the predictor stage


    !*******************************************************************************
    ! Step 2: Compute cell-centered reference and cell-corner electric fields  
    !*******************************************************************************    
    PRINT*,'Computing predictor CT cell center reference electric field.'

    CALL compute_cell_center_emf_2d(.FALSE.)  ! send false flag since this is the predictor stage

    PRINT*,'Computing predictor CT cell corner electric field.'

    CALL compute_cell_corner_emf_2d()


    !*******************************************************************************
    ! Step 3: Apply 1/2 dt CT update to cell-face magnetic fields
    !*******************************************************************************
    PRINT*,'Computing CT half update.'

    CALL update_CT(.FALSE.) ! send false flag since this is the predictor stage


    !*******************************************************************************
    ! Step 4: Apply dimensionally split 1/2 dt updates to density and momentum using
    !         transverse fluxes (including multidimensional source terms) immediately 
    !         followed by the "corrector" flux calculation.
    !*******************************************************************************
    PRINT*,'Computing state variables half update and final TVD fluxes in x direction.'

    ! 1/2 dt update along x
    CALL halfupdate_x_2d()
     
    ! compute corrector TVD x-fluxes     
    CALL compute_xfluxes_2d(.TRUE.)     ! send true flag since this is the corrector stage
        

    PRINT*,'Computing state variables half update and final TVD fluxes in y direction.'

    ! 1/2 dt update along y    
    CALL halfupdate_y_2d()
   
    ! compute corrector TVD y-fluxes     
    CALL compute_yfluxes_2d(.TRUE.)     ! send true flag since this is the corrector stage
    

    !*********************************************************************************
    ! Step 5: Compute the final cell-center reference and cell-corner electric fields   
    !*********************************************************************************
    PRINT*,'Computing corrector CT cell-center reference electric field.' 

    CALL compute_cell_center_emf_2d(.TRUE.)   ! send true flag since this is the corrector stage

    PRINT*,'Computing corrector CT cell-corner electric field.' 

    CALL compute_cell_corner_emf_2d()

    
    !*******************************************************************************
    ! Step 6: Apply full dt update to density, momentum and cell face magnetic 
    !         fields using the new fluxes and cell-corner electic fields
    !******************************************************************************* 
    PRINT*,'Computing final unsplit update of state variables and CT.'
    
    CALL update_state_vector_2d()
   
    ! copy original cell face magnetic fields values into main array
    bface_2d = bface
    
    CALL update_CT(.TRUE.) ! send true flag since this is the corrector stage

    PRINT*,'Update for time step completed.'

END SUBROUTINE solve_2d


! corrector stage half-dt update of state vector along x-direction using y-fluxes
SUBROUTINE halfupdate_x_2d()

   INTEGER :: i,j
   REAL*8 :: dtdx, idx, dBx, smomx, smomy, smomz, sBz

   dtdx = 0.5d0 * dt / dx   ! for half-step
   idx = 1.d0 / dx
   
   
   ! first, copy state vector into our work array
   q = q_2d
   
   ! now evolve the state vector by a half-stime step
   DO j = 1-nb+2, ny+nb-2
        DO i = 1-nb+2, nx+nb-2
   
            ! set the multi-dimensinal source terms
            dBx   = bface(i,j,1) - bface(i-1,j,1)   ! using initial cell-face magnetic field
            smomx = q(i,j,5) * dBx
            smomy = q(i,j,6) * dBx
            smomz = q(i,j,7) * dBx
            sBz   = (q(i,j,4) / q(i,j,1)) * dBx 
   
            ! evolve using the y-fluxes
            q(i,j,1) = q(i,j,1) - dtdx * ( Fs_buff(i,j,1,2) - Fs_buff(i,j-1,1,2) )
            q(i,j,2) = q(i,j,2) - dtdx * ( Fs_buff(i,j,2,2) - Fs_buff(i,j-1,2,2) ) + dtdx * smomx
            q(i,j,3) = q(i,j,3) - dtdx * ( Fs_buff(i,j,3,2) - Fs_buff(i,j-1,3,2) ) + dtdx * smomy
            q(i,j,4) = q(i,j,4) - dtdx * ( Fs_buff(i,j,4,2) - Fs_buff(i,j-1,4,2) ) + dtdx * smomz
            q(i,j,5) = q(i,j,5) - dtdx * ( Fs_buff(i,j,5,2) - Fs_buff(i,j-1,5,2) ) 
            !q(i,j,6) = q(i,j,6) - dtdx * ( Fs_buff(i,j,6,2) - Fs_buff(i,j-1,6,2) ) 
            q(i,j,7) = q(i,j,7) - dtdx * ( Fs_buff(i,j,7,2) - Fs_buff(i,j-1,7,2) ) + dtdx * sBz

        END DO
    END DO

END SUBROUTINE halfupdate_x_2d


! corrector stage half-dt update of state vector along y-direction using x-fluxes
SUBROUTINE halfupdate_y_2d()

   INTEGER :: i,j
   REAL*8 :: dtdx, idx, dBy, smomx, smomy, smomz, sBz

   dtdx = 0.5d0 * dt / dx   ! for half-step
   idx = 1.d0 / dx
   
   
    ! first rotate the stave vector array so that fastest index runs along y-direction
    ! (also need to swap x<->y components for momentum and B field)    
    DO i = 1-nb, ny+nb 
        DO j = 1-nb,nx+nb
        
            q(j,i,1) = q_2d(i,j,1)
            q(j,i,2) = q_2d(i,j,3)
            q(j,i,3) = q_2d(i,j,2)
            q(j,i,4) = q_2d(i,j,4)
            q(j,i,5) = q_2d(i,j,6)
            q(j,i,6) = q_2d(i,j,5)
            q(j,i,7) = q_2d(i,j,7)
         
        END DO
    END DO
   
    ! now evolve the state vector by a half-stime step
    ! (Remember,  q is already rotated, so need to be careful about that)
   DO i = 1-nb+2, ny+nb-2
        DO j = 1-nb+2, nx+nb-2
   
            ! set the multi-dimensinal source terms
            dBy   = bface(i,j,2) - bface(i,j-1,2)   ! using initial cell-face magnetic field
            smomx = q(j,i,5) * dBy
            smomy = q(j,i,6) * dBy
            smomz = q(j,i,7) * dBy
            sBz   = (q(j,i,4) / q(j,i,1)) * dBy 
   
            ! evolve using the use the x-fluxes
            ! (make sure to use the right flux components since the x,y components are swapped)
            q(j,i,1) = q(j,i,1) - dtdx * ( Fs_buff(i,j,1,1) - Fs_buff(i-1,j,1,1) )
            q(j,i,2) = q(j,i,2) - dtdx * ( Fs_buff(i,j,3,1) - Fs_buff(i-1,j,3,1) ) + dtdx * smomx
            q(j,i,3) = q(j,i,3) - dtdx * ( Fs_buff(i,j,2,1) - Fs_buff(i-1,j,2,1) ) + dtdx * smomy
            q(j,i,4) = q(j,i,4) - dtdx * ( Fs_buff(i,j,4,1) - Fs_buff(i-1,j,4,1) ) + dtdx * smomz
            q(j,i,5) = q(j,i,5) - dtdx * ( Fs_buff(i,j,6,1) - Fs_buff(i-1,j,6,1) ) 
            !q(j,i,6) = q(j,i,6) - dtdx * ( Fs_buff(i,j,5,1) - Fs_buff(i-1,j,5,1) ) 
            q(j,i,7) = q(j,i,7) - dtdx * ( Fs_buff(i,j,7,1) - Fs_buff(i-1,j,7,1) ) + dtdx * sBz

        END DO
    END DO

END SUBROUTINE halfupdate_y_2d


! Computes 1/2 dt (precictor) TVD x-fluxes
SUBROUTINE compute_xfluxes_2d(corrector_stage) 

    LOGICAL, INTENT(IN) :: corrector_stage
    INTEGER :: i,j
    REAL*8 :: dtdx


    dtdx = dt / dx   ! for half-step

    IF(.NOT. corrector_stage) THEN
        
        dtdx = 0.5d0 * dtdx ! for half-update during predictor stage

        ! first, copy state vector into our work array (don't do this if this is the corrector stage since it was already done during half update)
        q = q_2d
    
    END IF 
 
    PRINT*,'Computing Roe average state.' 
    
    ! now sweep through strips along x-direction and compute x-interface fluxes   
    CALL compute_roe_avg_state_1d(1, nx, 1, ny)
    
    PRINT*,'Computing eigenvalues.' 
    
    CALL compute_eigenvalues_1d(1, nx, 1, ny)
          
    PRINT*,'Computing eigenvectors.'
    
    CALL compute_eigenvectors_1d(1, nx, 1, ny)

    PRINT*,'Computing x-fluxes.'

    ! for predictor stage, store TVD fluxes in the buffer array
    IF(.NOT. corrector_stage) THEN
        
        ! clear TVD flux array
        Fs_buff = 0.d0
        
        CALL compute_TVD_fluxes_1d(Fs_buff, 1, nx, 1, ny, dtdx, 1)

        ! this is a good place to store the cell-face center electric fields
        DO j = 1-nb, ny+nb
            DO i = 1-nb+1, nx+nb-2    
                emf_2d(i,j,1) = -Fs_buff(i,j,6,1)   ! this is Ez_i+1/2,j (i.e. negative By component of the x-flux)
            END DO
        END DO    
    
    ELSE
    
        ! clear TVD flux array
        Fs = 0.d0
        
        ! for corrector stage, we store fluxes in the main array
        CALL compute_TVD_fluxes_1d(Fs, 1, nx, 1, ny, dtdx, 1)

        ! this is a good place to store the cell-face center electric fields
        DO j = 1-nb, ny+nb
            DO i = 1-nb+1, nx+nb-2   
                emf_2d(i,j,1) = -Fs(i,j,6,1)   ! this is Ez_i+1/2,j (i.e. negative By component of the x-flux)
            END DO
        END DO    
    
    END IF
    

END SUBROUTINE compute_xfluxes_2d 


! Computes 1/2 dt (predictor) TVD y-fluxes
SUBROUTINE compute_yfluxes_2d(corrector_stage) 

    LOGICAL, INTENT(IN) :: corrector_stage
    INTEGER :: i,j
    REAL*8 :: dtdx


    dtdx = dt / dx      
    
    IF(.NOT. corrector_stage) THEN
    
        dtdx = 0.5d0 * dtdx ! for half-update during predictor stage
    
        ! first rotate the stave vector array so that fastest index runs along y-direction (don't do this for corrector stage since it was already done during half update)
        ! (also need to swap x<->y components for momentum and B field)
        DO i = 1-nb, ny+nb 
            DO j = 1-nb, nx+nb
        
                q(j,i,1) = q_2d(i,j,1)
                q(j,i,3) = q_2d(i,j,2)
                q(j,i,2) = q_2d(i,j,3)
                q(j,i,4) = q_2d(i,j,4)
                q(j,i,6) = q_2d(i,j,5)
                q(j,i,5) = q_2d(i,j,6)
                q(j,i,7) = q_2d(i,j,7)
         
            END DO
        END DO
     
     END IF
   
     PRINT*,'Computing Roe average state.' 

     ! now sweep through strips along y-direction and compute y-interface fluxes   
     CALL compute_roe_avg_state_1d(1, ny, 1, nx)
    
     PRINT*,'Computing eigenvalues.'
     
     CALL compute_eigenvalues_1d(1, ny, 1, nx)
       
     PRINT*,'Computing eigenvectors.'
    
     CALL compute_eigenvectors_1d(1, ny, 1, nx)

     PRINT*,'Computing y-fluxes.'

    ! for predictor stage, store TVD fluxes in the buffer array
    IF(.NOT. corrector_stage) THEN
    
        ! clear TVD flux array
        Fs = 0.d0 
        
        CALL compute_TVD_fluxes_1d(Fs, 1, ny, 1, nx, dtdx, 2)


        ! Rotate the y-fluxes in flux buffer and store them in main flux array 
        DO j = 1-nb, ny+nb 
            DO i = 1-nb, nx+nb
        
                Fs_buff(i,j,1,2) = Fs(j,i,1,2)
                Fs_buff(i,j,3,2) = Fs(j,i,2,2)
                Fs_buff(i,j,2,2) = Fs(j,i,3,2)
                Fs_buff(i,j,4,2) = Fs(j,i,4,2)
                Fs_buff(i,j,6,2) = Fs(j,i,5,2)
                Fs_buff(i,j,5,2) = Fs(j,i,6,2)         
                Fs_buff(i,j,7,2) = Fs(j,i,7,2)         
            
            END DO
        END DO
    
        ! this is a good place to store the cell-face center electric fields
        DO j = 1-nb+1, ny+nb-2
            DO i = 1-nb, nx+nb   
                emf_2d(i,j,2) = Fs_buff(i,j,5,2)   ! this is Ez_i,j+1/2 (Note: After rotation, component 5 of the y-flux is Ez_i,j+1/2)
            END DO
        END DO      

    ELSE 
    
        ! clear TVD flux array
        Fs_buff = 0.d0
    
        CALL compute_TVD_fluxes_1d(Fs_buff, 1, ny, 1, nx, dtdx, 2)


        ! Rotate the y-fluxes in flux buffer and store them in main flux array 
        DO j = 1-nb, ny+nb 
            DO i = 1-nb,nx+nb
        
                Fs(i,j,1,2) = Fs_buff(j,i,1,2)
                Fs(i,j,3,2) = Fs_buff(j,i,2,2)
                Fs(i,j,2,2) = Fs_buff(j,i,3,2)
                Fs(i,j,4,2) = Fs_buff(j,i,4,2)
                Fs(i,j,6,2) = Fs_buff(j,i,5,2)
                Fs(i,j,5,2) = Fs_buff(j,i,6,2)         
                Fs(i,j,7,2) = Fs_buff(j,i,7,2)         
            
            END DO
        END DO
    
        ! this is a good place to store the cell-face center electric fields
        DO j = 1-nb+1, ny+nb-2
            DO i = 1-nb, nx+nb  
                emf_2d(i,j,2) = Fs(i,j,5,2)   ! this is Ez_i,j+1/2 (Note: After rotation, component 5 of the y-flux is Ez_i,j+1/2)
            END DO
        END DO      

    
    END IF


END SUBROUTINE compute_yfluxes_2d 


! Compute the cell-centered refernece emf: Ez_i,j

! For the corrector stage, need to compute vx, vy at half-time step by doing 
! an auxiliary half-update of the density and momentum x,y variables. For this half-update
! must use the predictor fluxes, not final corrector fluxes. NEED TO FIX THIS, BECAUSE
! CURRENTLY WE'RE USING THE WRONG FLUXES (even though it works) !!! 

SUBROUTINE compute_cell_center_emf_2d(corrector_stage)

    LOGICAL, INTENT(IN) :: corrector_stage
    INTEGER :: i, j
    REAL*8 :: rho, vx, vy, momx, momy, Bx, By, dtdx
    
    
    ! For the corrector stage, wee need the half-updated cell-center momentum
    IF(corrector_stage) THEN
            
        dtdx = 0.5d0 * dt / dx
        
        DO j = 1-nb+2, ny+nb-2
            DO i = 1-nb+2, nx+nb-2     

                ! To get cell-centered emf, we need cell-centered velocity at half-time step.
                ! Can get those by applying an unsplit conservative half-update of density and 
                ! momentum using the already available 1/2 dt fluxes. 
                rho  = q_2d(i,j,1) - dtdx * (Fs(i,j,1,1) - Fs(i-1,j,1,1) + Fs(i,j,1,2) - Fs(i,j-1,1,2))  
                momx = q_2d(i,j,2) - dtdx * (Fs(i,j,2,1) - Fs(i-1,j,2,1) + Fs(i,j,2,2) - Fs(i,j-1,2,2))  
                momy = q_2d(i,j,3) - dtdx * (Fs(i,j,3,1) - Fs(i-1,j,3,1) + Fs(i,j,3,2) - Fs(i,j-1,3,2))  
        
                vx = momx / rho
                vy = momy / rho
         
                ! Get magnetic field components at half-time step by averaging over cell-face
                ! values obtained from step 3.  
                Bx = 0.5d0 * (bface_2d(i,j,1) + bface_2d(i-1,j,1))         
                By = 0.5d0 * (bface_2d(i,j,2) + bface_2d(i,j-1,2))    
        
                emf_2d(i,j,3) = Bx * vy - By * vx          
        
            END DO
        END DO    
        
    ELSE
        
        DO j = 1-nb+2, ny+nb-2
            DO i = 1-nb+2, nx+nb-2     

                emf_2d(i,j,3) = ( q_2d(i,j,5) * q_2d(i,j,3)  - q_2d(i,j,6) * q_2d(i,j,2) ) / q_2d(i,j,1)       
            
            END DO
        END DO    
        
    END IF

END SUBROUTINE compute_cell_center_emf_2d


! Compute the cell-corner emf
SUBROUTINE compute_cell_corner_emf_2d()

    INTEGER :: i, j
    REAL*8 :: vx, vy, idx2, &
              dyL_11_14, dyR_11_14, dyL_11_34, dyR_11_34, &
              dxL_14_11, dxR_14_11, dxL_34_11, dxR_34_11, &
              dy_12_14, dy_12_34, dx_14_12, dx_34_12
    
    idx2 = 2.d0 / dx
    
    DO j = 1-nb+2, ny+nb-3
        DO i = 1-nb+2, nx+nb-3  
        
            ! get normal velocity at x-interface
            vx = 0.5d0 * (q_2d(i,j,2)/q_2d(i,j,1) + q_2d(i+1,j,2)/q_2d(i+1,j,1) )
            
            ! get normal velocity at y-interface
            vy = 0.5d0 * (q_2d(i,j,3)/q_2d(i,j,1) + q_2d(i,j+1,3)/q_2d(i,j+1,1) )
            
            ! compute upwinded interface gradients of electric field at x-interface
            dyL_11_14 = idx2 * (emf_2d(i,j,2) - emf_2d(i,j,3)) 
            dyR_11_14 = idx2 * (emf_2d(i+1,j,2) - emf_2d(i+1,j,3)) 
            dyL_11_34 = idx2 * (emf_2d(i,j+1,3) - emf_2d(i,j,2))
            dyR_11_34 = idx2 * (emf_2d(i+1,j+1,3) - emf_2d(i+1,j,2))

            IF(vx .GT. ZERO) THEN         
            
                dy_12_14 =  dyL_11_14
                dy_12_34 =  dyL_11_34 
                
            ELSE IF(vx .LT. -ZERO) THEN
            
                dy_12_14 =  dyR_11_14
                dy_12_34 =  dyR_11_34 

            ELSE
            
                dy_12_14 = 0.5d0 * (dyL_11_14 + dyR_11_14)
                dy_12_34 = 0.5d0 * (dyL_11_34 + dyR_11_34)
                
            END IF
            
            ! compute upwinded interface gradients of electric field at y-interface
            dxL_14_11 = idx2 * (emf_2d(i,j,1) - emf_2d(i,j,3))             
            dxR_14_11 = idx2 * (emf_2d(i,j+1,1) - emf_2d(i,j+1,3))             
            dxL_34_11 = idx2 * (emf_2d(i+1,j,3) - emf_2d(i,j,1))
            dxR_34_11 = idx2 * (emf_2d(i+1,j+1,3) - emf_2d(i,j+1,1))
                        

            IF(vy .GT. ZERO) THEN         
            
                dx_14_12 =  dxL_14_11
                dx_34_12 =  dxL_34_11 
                
            ELSE IF(vy .LT. -ZERO) THEN
            
                dx_14_12 =  dxR_14_11
                dx_34_12 =  dxR_34_11 

            ELSE
            
                dx_14_12 = 0.5d0 * (dxL_14_11 + dxR_14_11)
                dx_34_12 = 0.5d0 * (dxL_34_11 + dxR_34_11)
                
            END IF
            
        
            ! cell-corner electric field z component: Ez_i+1/2,j+1/2 
            emf_corner(i,j) =  0.25d0 * (emf_2d(i,j,1) + emf_2d(i,j,2) + emf_2d(i,j+1,1) + emf_2d(i+1,j,2)) + &
                               0.125d0 * dx * ((dy_12_14 - dy_12_34) + (dx_14_12 - dx_34_12)) 
                        
        END DO
    END DO    
        

END SUBROUTINE compute_cell_corner_emf_2d


! updates cell-faces magnetic fields via CT
SUBROUTINE update_CT(corrector_stage)

    LOGICAL, INTENT(IN) :: corrector_stage
    REAL*8 :: dtdx
    INTEGER :: i, j

    dtdx = dt / dx

    ! for the initial predictor stage, we only apply a half-time update
    IF(.NOT. corrector_stage) dtdx = 0.5d0 * dtdx
    
    
    DO j = 1-nb+3, ny+nb-3
        DO i = 1-nb+2, nx+nb-3  

            ! Bx_i+1/2,j
            bface_2d(i,j,1) = bface_2d(i,j,1) - dtdx * (emf_corner(i,j) - emf_corner(i,j-1))  
                
               
        END DO
    END DO


    DO j = 1-nb+2, ny+nb-3
        DO i = 1-nb+3, nx+nb-3  
                
            ! By_i,j+1/2
            bface_2d(i,j,2) = bface_2d(i,j,2) + dtdx * (emf_corner(i,j) - emf_corner(i-1,j))                
               
        END DO
    END DO


    ! If this is the final (i.e. corrector stage) CT update, then also compute the cell-centered magnetic fields
    ! by averaging over the cell-face values (i.e. linear interpolation)
    IF(corrector_stage) THEN
    
        DO j = 1-nb+3, ny+nb-3
            DO i = 1-nb+3, nx+nb-3  

                ! Bx_i,j
                q_2d(i,j,5) = 0.5d0 * (bface_2d(i,j,1) + bface_2d(i-1,j,1))  
                
                ! By_i,j
                q_2d(i,j,6) = 0.5d0 * (bface_2d(i,j,2) + bface_2d(i,j-1,2))                
               
            END DO
        END DO
 
 
    END IF

END SUBROUTINE update_CT



! Compute cell-interface roe averaged state
SUBROUTINE compute_roe_avg_state_1d(ilow, ihi, jlow, jhi)

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi
    REAL*8 :: rho2, vx2, vy2, vz2, Bx2, By2, Bz2, isrho
    INTEGER :: i, j
    
       
    DO j = 1-nb, jhi+nb     
    DO i = 1-nb, ihi+nb-1

        ! compute Roe-averaged primitive variables (this is just an arithmetic average)
        rho2 = 0.5d0 * ( q(i,j,1) + q(i+1,j,1) )
        isrho = 1.D0 / SQRT(rho2)
        vx2  = 0.5d0 * ( q(i,j,2)/q(i,j,1) + q(i+1,j,2)/q(i+1,j,1) )
        vy2  = 0.5d0 * ( q(i,j,3)/q(i,j,1) + q(i+1,j,3)/q(i+1,j,1) )
        vz2  = 0.5d0 * ( q(i,j,4)/q(i,j,1) + q(i+1,j,4)/q(i+1,j,1) )
        Bx2  = 0.5d0 * ( q(i,j,5) + q(i+1,j,5) )
        By2  = 0.5d0 * ( q(i,j,6) + q(i+1,j,6) )
        Bz2  = 0.5d0 * ( q(i,j,7) + q(i+1,j,7) )
        
        qavg(i,j,1) = rho2
        qavg(i,j,2) = vx2
        qavg(i,j,3) = vy2
        qavg(i,j,4) = vz2
        qavg(i,j,5) = Bx2 * isrho
        qavg(i,j,6) = By2 * isrho
        qavg(i,j,7) = Bz2 * isrho           

    END DO 
    END DO


END SUBROUTINE compute_roe_avg_state_1d


! Compute eigenvalues
SUBROUTINE compute_eigenvalues_1d(ilow, ihi, jlow, jhi)

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi
    REAL*8 :: rho2, vx2, vy2, vz2, bx2, by2, bz2, asqr
    REAL*8 :: cf, ca, cs, bsqr, btsqr, bx, tmp1, tmp2, isrho
    INTEGER :: i, j
        
    DO j = 1-nb, jhi+nb     
    DO i = 1-nb, ihi+nb-1

        rho2  = qavg(i,j,1)
        vx2   = qavg(i,j,2)
        asqr  = sound_speed**2
        bx    = qavg(i,j,5)
        bsqr  = qavg(i,j,5)**2 + qavg(i,j,6)**2 + qavg(i,j,7)**2
        !btsqr  =  qavg(i,j,6)**2 + qavg(i,j,7)**2
        tmp1  = asqr + bsqr        
        !tmp2  = SQRT( (asqr - bsqr)**2 + 4.d0 * asqr * btsqr )
        tmp2  = SQRT( (asqr + bsqr)**2 - 4.d0 * asqr * bx**2 )

        cs = SQRT( MAX(0.D0, 0.5d0 * (tmp1 - tmp2) ))  ! slow mode speed
        ca = ABS(bx)                                   ! alfven speed
        cf = SQRT( MAX(0.D0, 0.5d0 * (tmp1 + tmp2) ))  ! fast mode speed
        
        ! wave speeds
        wave_speeds(i,j,1) = cs
        wave_speeds(i,j,2) = ca
        wave_speeds(i,j,3) = cf
        
        ! eigenvalues
        eigenvalues(i,j,1) = vx2 - cf
        eigenvalues(i,j,2) = vx2 - ca
        eigenvalues(i,j,3) = vx2 - cs
        eigenvalues(i,j,4) = vx2 + cs
        eigenvalues(i,j,5) = vx2 + ca
        eigenvalues(i,j,6) = vx2 + cf


        !IF(cf .LT. ca) THEN
        !    PRINT*,'Round off error. cf - ca = ',cf-ca
        !    STOP        
        !END IF

    END DO 
    END DO
    
    
END SUBROUTINE compute_eigenvalues_1d


! Compute eigenvectors
SUBROUTINE compute_eigenvectors_1d(ilow, ihi, jlow, jhi)

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi
    REAL*8 :: rho2, vx2, vy2, vz2, bx2, by2, bz2, a1, asqr, isrho, amc
    REAL*8 :: cf, ca, cs, bsqr, bt, tmp1, tmp2, tmp3, tmp4, tmp5, signbx, signbt, itheta1, itheta2
    REAL*8 :: alphas(1-nb:nmax+nb,1-nb:nmax+nb), alphaf(1-nb:nmax+nb,1-nb:nmax+nb), &
              betay(1-nb:nmax+nb,1-nb:nmax+nb), betaz(1-nb:nmax+nb,1-nb:nmax+nb)
    REAL*8 :: delq(6), Lminus(6), Lplus(6)
    INTEGER :: i, j
     
     
    !PRINT*,'*** CHECKPOINT 1 ***' 
         
    DO j = 1-nb, jhi+nb     
    DO i = 1-nb, ihi+nb-1

        asqr   = sound_speed**2      
        bt = SQRT(qavg(i,j,6)**2 + qavg(i,j,7)**2) ! magnitude of B-transverse

        ! compute eigenvector renormalization co-efficients, taking limiting values for degenracies (K98 eqn. 2.13-15)
        IF(ABS(wave_speeds(i,j,3)**2 - wave_speeds(i,j,1)**2) .LT. ZERO) THEN
            alphaf(i,j) = 1.D0
            alphas(i,j) = 1.D0
        ELSE
            tmp1 = 1.D0 / SQRT( MAX(0.d0, wave_speeds(i,j,3)**2 - wave_speeds(i,j,1)**2) )
            alphaf(i,j) = SQRT( MAX(0.d0, wave_speeds(i,j,3)**2 - wave_speeds(i,j,2)**2) ) * tmp1      ! Need to be careful with square roots. Sometimes round-off errors can make cf < ca, or ca < cs.   
            alphas(i,j) = SQRT( MAX(0.d0, wave_speeds(i,j,3)**2 - asqr) ) * tmp1
        END IF
 
        IF(bt .LT. ZERO) THEN        
            betay(i,j) = 1.D0 / SQRT(2.D0)
            betaz(i,j) = 1.D0 / SQRT(2.D0)  
        ELSE
            tmp2 = 1.D0 / bt
            betay(i,j) = qavg(i,j,6) * tmp2
            betaz(i,j) = qavg(i,j,7) * tmp2       
        END IF
        
    END DO    
    END DO


    !PRINT*,'*** CHECKPOINT 2 ***'   
    
    !*************************************************
    ! compute the right eigenvectors (K98 eq 2.12)
    !*************************************************
    
    DO j = 1-nb, jhi+nb     
    DO i = 1-nb, ihi+nb-1
    
        isrho = 1.D0 / SQRT(qavg(i,j,1))       
        bt = SQRT(qavg(i,j,6)**2 + qavg(i,j,7)**2) ! magnitude of B-transverse
        amc = sound_speed**2 - wave_speeds(i,j,2)**2
        tmp1 = alphaf(i,j) * qavg(i,j,3)
        tmp2 = alphas(i,j) * betay(i,j) * qavg(i,j,5)
        tmp3 = alphaf(i,j) * qavg(i,j,4)
        tmp4 = alphas(i,j) * betaz(i,j) * qavg(i,j,5)
        tmp5 = alphas(i,j) * wave_speeds(i,j,3) * isrho
        
        ! fast mode minus eigenvector 
        Rk(i,j,1,1) = alphaf(i,j)
        Rk(i,j,2,1) = alphaf(i,j) * eigenvalues(i,j,1) 
        Rk(i,j,3,1) = tmp1 + tmp2
        Rk(i,j,4,1) = tmp3 + tmp4 
        Rk(i,j,5,1) = tmp5 * betay(i,j)
        Rk(i,j,6,1) = tmp5 * betaz(i,j)
        
        ! fast mode plus eigenvector 
        Rk(i,j,1,6) = alphaf(i,j)
        Rk(i,j,2,6) = alphaf(i,j) * eigenvalues(i,j,6) 
        Rk(i,j,3,6) = tmp1 - tmp2
        Rk(i,j,4,6) = tmp3 - tmp4 
        Rk(i,j,5,6) = tmp5 * betay(i,j)
        Rk(i,j,6,6) = tmp5 * betaz(i,j)
        

        ! check for continuity (K98 eqn. 2.20)
        IF(amc .LT. ZERO) THEN
        
            IF(ABS(qavg(i,j,6)) .GE. ZERO) THEN
                signbt = SIGN(1.D0, qavg(i,j,6))
            ELSE
                signbt = SIGN(1.D0, qavg(i,j,7))            
            END IF
            
            Rk(i,j,:,1) = signbt * Rk(i,j,:,1) 
            Rk(i,j,:,6) = signbt * Rk(i,j,:,6)
        
        END IF

    END DO
    END DO 

    DO j = 1-nb, jhi+nb     
    DO i = 1-nb, ihi+nb-1
    
        isrho = 1.D0 / SQRT(qavg(i,j,1))       
        signbx = SIGN(1.D0, qavg(i,j,5))
        
        ! alfven mode minus eigenvector 
        Rk(i,j,1,2) = 0.D0
        Rk(i,j,2,2) = 0.D0 
        Rk(i,j,3,2) = signbx * betaz(i,j)
        Rk(i,j,4,2) = -signbx * betay(i,j)
        Rk(i,j,5,2) = isrho * betaz(i,j)
        Rk(i,j,6,2) = -isrho * betay(i,j)
        
        ! alfven mode plus eigenvector 
        Rk(i,j,1,5) = 0.D0
        Rk(i,j,2,5) = 0.D0 
        Rk(i,j,3,5) = -signbx * betaz(i,j)
        Rk(i,j,4,5) = signbx * betay(i,j)
        Rk(i,j,5,5) = isrho * betaz(i,j)
        Rk(i,j,6,5) = -isrho * betay(i,j)
        
    END DO    
    END DO 
    
    DO j = 1-nb, jhi+nb     
    DO i = 1-nb, ihi+nb-1
    
        isrho = 1.D0 / SQRT(qavg(i,j,1))       
        signbx = SIGN(1.D0, qavg(i,j,5))
        amc = sound_speed**2 - wave_speeds(i,j,2)**2
        tmp1 = alphas(i,j) * qavg(i,j,3)
        tmp2 = alphaf(i,j) * betay(i,j) * sound_speed * signbx 
        tmp3 = alphas(i,j) * qavg(i,j,4)
        tmp4 = alphaf(i,j) * betaz(i,j) * sound_speed * signbx
        tmp5 = alphaf(i,j) * (sound_speed**2) * isrho /  wave_speeds(i,j,3)
        
        ! slow mode minus eigenvector 
        Rk(i,j,1,3) = alphas(i,j)
        Rk(i,j,2,3) = alphas(i,j) * eigenvalues(i,j,3) 
        Rk(i,j,3,3) = tmp1 - tmp2
        Rk(i,j,4,3) = tmp3 - tmp4 
        Rk(i,j,5,3) = -tmp5 * betay(i,j)
        Rk(i,j,6,3) = -tmp5 * betaz(i,j)
        
        ! slow mode plus eigenvector 
        Rk(i,j,1,4) = alphas(i,j)
        Rk(i,j,2,4) = alphas(i,j) * eigenvalues(i,j,4) 
        Rk(i,j,3,4) = tmp1 + tmp2
        Rk(i,j,4,4) = tmp3 + tmp4 
        Rk(i,j,5,4) = -tmp5 * betay(i,j)
        Rk(i,j,6,4) = -tmp5 * betaz(i,j)
        
        ! check for continuity (K98 eqn. 2.20)
        IF(amc .GT. ZERO) THEN
        
            IF(ABS(qavg(i,j,6)) .GE. ZERO) THEN
                signbt = SIGN(1.D0, qavg(i,j,6))
            ELSE
                signbt = SIGN(1.D0, qavg(i,j,7))            
            END IF
        
            Rk(i,j,:,3) = signbt * Rk(i,j,:,3) 
            Rk(i,j,:,4) = signbt * Rk(i,j,:,4)
        
        END IF
       
    END DO     
    END DO
    

    !PRINT*,'*** CHECKPOINT 3 ***' 


    !*****************************************************
    ! compute the characteristic variables (K98 eqn. 2.27)
    !*****************************************************
    
    DO j = 1-nb, jhi+nb     
    DO i = 1-nb, ihi+nb-1
                            
        signbx = SIGN(1.D0, qavg(i,j,5))
        amc = sound_speed**2 - wave_speeds(i,j,2)**2

        ! qavg_i+1 - qavg_i (Note: The Bx component has to be excluded)
        delq(1:4) = q(i+1,j,1:4) - q(i,j,1:4) 
        delq(5) = q(i+1,j,6) - q(i,j,6)
        delq(6) = q(i+1,j,7) - q(i,j,7)
 
        itheta1 = 1.D0 / (2.D0 * ( alphaf(i,j)**2 * sound_speed**2 + alphas(i,j)**2 * wave_speeds(i,j,3)**2))  
        itheta2 = 1.D0 / (2.D0 * ( alphaf(i,j)**2 * sound_speed * wave_speeds(i,j,3) + &
                  alphas(i,j)**2 * wave_speeds(i,j,1) * wave_speeds(i,j,2) ) )
       
        tmp1 = itheta1 * alphaf(i,j) * sound_speed**2
        tmp2 = itheta2 * ( - alphaf(i,j) * sound_speed * qavg(i,j,2) + &
               alphas(i,j) * wave_speeds(i,j,1) * signbx * ( betay(i,j) * qavg(i,j,3) + betaz(i,j) * qavg(i,j,4)) )         
 
        tmp3 = itheta2 * alphas(i,j) * wave_speeds(i,j,1) * signbx 
        tmp4 = itheta1 * alphas(i,j) * wave_speeds(i,j,3) * SQRT(qavg(i,j,1))    
        
        ! Fast mode minus and plus eigenvector components (K98 eqn. 2.16-19)
        Lminus(1) =  tmp1 - tmp2
        Lplus(1)  = tmp1 + tmp2
        
        Lminus(2) = -itheta2 * alphaf(i,j) * sound_speed  
        Lplus(2)  = -Lminus(2)
       
        Lminus(3) = tmp3 * betay(i,j)   
        Lplus(3)  = -Lminus(3)
       
        Lminus(4) = tmp3 * betaz(i,j)
        Lplus(4)  = -Lminus(4)
       
        Lminus(5) = tmp4 * betay(i,j)
        Lplus(5)  = Lminus(5)     
        
        Lminus(6) = tmp4 * betaz(i,j)
        Lplus(6)  = Lminus(6)     
        
         ! check for continuity (K98 eqn. 2.20)
        IF(amc .LT. ZERO) THEN
        
            IF(ABS(qavg(i,j,6)) .GE. ZERO) THEN
                signbt = SIGN(1.D0, qavg(i,j,6))
            ELSE
                signbt = SIGN(1.D0, qavg(i,j,7))            
            END IF
            
            Lminus = signbt * Lminus
            Lplus  = signbt * Lplus
        
        END IF
        
        
        ! Fast mode minus and plus chacteristic variables
        ck(i,j,1) = delq(1) * Lminus(1) + delq(2) * Lminus(2) + delq(3) * Lminus(3) + & 
                    delq(4) * Lminus(4) + delq(5) * Lminus(5) + delq(6) * Lminus(6)  

        ck(i,j,6) = delq(1) * Lplus(1) + delq(2) * Lplus(2) + delq(3) * Lplus(3) + & 
                    delq(4) * Lplus(4) + delq(5) * Lplus(5) + delq(6) * Lplus(6)
        
    END DO 
    END DO
    
    DO j = 1-nb, jhi+nb     
    DO i = 1-nb, ihi+nb-1
              
        signbx = SIGN(1.D0, qavg(i,j,5))

        ! qavg_i+1 - qavg_i (Note: The Bx component has to be excluded)
        delq(1:4) = q(i+1,j,1:4) - q(i,j,1:4) 
        delq(5) = q(i+1,j,6) - q(i,j,6)
        delq(6) = q(i+1,j,7) - q(i,j,7)
        
        tmp1 = 0.5d0 * (betaz(i,j) * qavg(i,j,3) - betay(i,j) * qavg(i,j,4)) * signbx
        tmp2 = 0.5d0 * signbx
        tmp3 = 0.5d0 * SQRT(qavg(i,j,1))
        
        ! Alfven mode minus and plus eigenvector components (K98 eqn. 2.16-19)
        Lminus(1) =  -tmp1
        Lplus(1)  = tmp1 
        
        Lminus(2) = 0.D0  
        Lplus(2)  = 0.D0
       
        Lminus(3) = tmp2 * betaz(i,j)
        Lplus(3)  = -Lminus(3)
       
        Lminus(4) = -tmp2 * betay(i,j)
        Lplus(4)  = -Lminus(4)
       
        Lminus(5) = tmp3 * betaz(i,j)
        Lplus(5)  = Lminus(5)     
        
        Lminus(6) = -tmp3 * betay(i,j)
        Lplus(6)  = Lminus(6)     
        
        
        ! Alfven mode minus and plus chacteristic variables       
        ck(i,j,2) = delq(1) * Lminus(1) + delq(2) * Lminus(2) + delq(3) * Lminus(3) + & 
                    delq(4) * Lminus(4) + delq(5) * Lminus(5) + delq(6) * Lminus(6)  

        ck(i,j,5) = delq(1) * Lplus(1) + delq(2) * Lplus(2) + delq(3) * Lplus(3) + & 
                    delq(4) * Lplus(4) + delq(5) * Lplus(5) + delq(6) * Lplus(6)
                  
    END DO 
    END DO
    
        
    DO j = 1-nb, jhi+nb     
    DO i = 1-nb, ihi+nb-1
              
        signbx = SIGN(1.D0, qavg(i,j,5))
        amc = sound_speed**2 - wave_speeds(i,j,2)**2 

        ! qavg_i+1 - qavg_i (Note: The Bx component has to be excluded)
        delq(1:4) = q(i+1,j,1:4) - q(i,j,1:4) 
        delq(5) = q(i+1,j,6) - q(i,j,6)
        delq(6) = q(i+1,j,7) - q(i,j,7)
        
        itheta1 = 1.D0 / (2.D0 * ( alphaf(i,j)**2 * sound_speed**2 + alphas(i,j)**2 * wave_speeds(i,j,3)**2))  
        itheta2 = 1.D0 / (2.D0 * ( alphaf(i,j)**2 * sound_speed * wave_speeds(i,j,3) + &
                  alphas(i,j)**2 * wave_speeds(i,j,1) * wave_speeds(i,j,2) ) )
       
        tmp1 = itheta1 * alphas(i,j) * wave_speeds(i,j,3)**2
        tmp2 = itheta2 * ( alphas(i,j) * wave_speeds(i,j,2) * qavg(i,j,2) + &
               alphaf(i,j) * wave_speeds(i,j,3) * signbx * ( betay(i,j) * qavg(i,j,3) + betaz(i,j) * qavg(i,j,4)) )         
 
        tmp3 = itheta2 * alphaf(i,j) * wave_speeds(i,j,3) * signbx 
        tmp4 = itheta1 * alphaf(i,j) * wave_speeds(i,j,3) * SQRT(qavg(i,j,1))    
        
        ! Slow mode minus and plus eigenvector components (K98 eqn. 2.16-19)
        Lminus(1) = tmp1 + tmp2
        Lplus(1)  = tmp1 - tmp2
        
        Lminus(2) = -itheta2 * alphas(i,j) * wave_speeds(i,j,2)  
        Lplus(2)  = -Lminus(2)
       
        Lminus(3) = -tmp3 * betay(i,j)   
        Lplus(3)  = -Lminus(3)
       
        Lminus(4) = -tmp3 * betaz(i,j)
        Lplus(4)  = -Lminus(4)
       
        Lminus(5) = -tmp4 * betay(i,j)
        Lplus(5)  = Lminus(5)     
        
        Lminus(6) = -tmp4 * betaz(i,j)
        Lplus(6)  = Lminus(6)     
        
         ! check for continuity (K98 eqn. 2.20)
        IF(amc .GT. ZERO) THEN
        
            IF(ABS(qavg(i,j,6)) .GE. ZERO) THEN
                signbt = SIGN(1.D0, qavg(i,j,6))
            ELSE
                signbt = SIGN(1.D0, qavg(i,j,7))            
            END IF
            
            Lminus = signbt * Lminus
            Lplus  = signbt * Lplus
        
        END IF
        
        
        ! Slow mode minus and plus chacteristic variables
        ck(i,j,3) = delq(1) * Lminus(1) + delq(2) * Lminus(2) + delq(3) * Lminus(3) + & 
                    delq(4) * Lminus(4) + delq(5) * Lminus(5) + delq(6) * Lminus(6)  

        ck(i,j,4) = delq(1) * Lplus(1) + delq(2) * Lplus(2) + delq(3) * Lplus(3) + & 
                    delq(4) * Lplus(4) + delq(5) * Lplus(5) + delq(6) * Lplus(6)
                  
    END DO 
    END DO
    

    !PRINT*,'*** CHECKPOINT 4 ***' 
    
    
END SUBROUTINE compute_eigenvectors_1d


! Compute the time-averaged cell-interface fluxes
SUBROUTINE compute_TVD_fluxes_1d(Flux,ilow, ihi, jlow, jhi, dtdx, fcomp)

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi, fcomp
    REAL*8, INTENT(IN) :: dtdx
    REAL*8, INTENT(INOUT) :: Flux(1-nb:nmax+nb,1-nb:nmax+nb,7,2)
    INTEGER :: i, j, k
    REAL*8 :: dxdt, gammak, betak
    REAL*8 :: Qk, chi
    REAL*8 :: gtilde1, gtilde2, Qk1, Qk2, chi1, chi2, signg
    REAL*8 :: gk(1-nb:nmax+nb, 1-nb:nmax+nb, nwaves), F(1-nb:ihi+nb, 1-nb:jhi+nb,nwaves)
    REAL*8 :: rho, vx, vy, vz, Bx, By, Bz


    dxdt = dx/dt
    
    ! compute cell center fluxes (instead of storing this in an array, maybe compute these ad hoc inside the TVD flux calculation loop?)
    DO j = 1-nb, jhi+nb     
    DO i = 1-nb, ihi+nb
    
        rho = q(i,j,1)
        vx = q(i,j,2) / rho
        vy = q(i,j,3) / rho
        vz = q(i,j,4) / rho
        Bx = q(i,j,5)
        By = q(i,j,6)
        Bz = q(i,j,7)
    
        F(i,j,1) = rho * vx
        F(i,j,2) = rho * (vx**2 + sound_speed**2) + 0.5d0 * (By**2 + Bz**2 - Bx**2)
        F(i,j,3) = rho * vx * vy - Bx * By
        F(i,j,4) = rho * vx * vz - Bx * Bz
        F(i,j,5) = vx * By - vy * Bx
        F(i,j,6) = vx * Bz - vz * Bx   
    
    END DO
    END DO

    ! compute the flux limiters
     DO k = 1, nwaves
        DO j = 1-nb, jhi+nb     
        DO i = 1-nb+1, ihi+nb-1
                
            chi1 = dtdx*eigenvalues(i-1,j,k) 
            chi2 = dtdx*eigenvalues(i,j,k) 
        
            ! Qk_i-1/2 (K98 eqn. 2.31)
            IF(ABS(chi1) .LT. (2.D0*eps(k))) THEN
                Qk1 = 0.25D0 * (chi1**2 / eps(k)) + eps(k)
            ELSE
                Qk1 = ABS(chi1)
            END IF    

            ! Qk_i+1/2 (K98 eqn. 2.31)
            IF(ABS(chi2) .LT. (2.D0*eps(k))) THEN
                Qk2 = 0.25D0 * (chi2**2 / eps(k)) + eps(k)
            ELSE
                Qk2 = ABS(chi2)
            END IF  
            
            !gk_tilde_i-1/2 (K98 eqn. 2.30)
            gtilde1 = 0.5d0 * (Qk1 - chi1**2) * ck(i-1,j,k) 
            
            !gk_tilde_i+1/2 (K98 eqn. 2.30)
            gtilde2 = 0.5d0 * (Qk2 - chi2**2) * ck(i,j,k) 
         
            signg = SIGN(1.D0, gtilde2)             
            
            ! flux limiter (minmod limiter)
            gk(i,j,k) = signg * MAX( 0.d0, MIN(ABS(gtilde2), gtilde1*signg) )         
            
            !gk(i,k) = signg * MAX( 0.d0, MIN(ABS(gtilde2), limiter_beta*gtilde1*signg), MIN(limiter_beta*ABS(gtilde2), gtilde1*signg) )   ! Sweby limiter      
                
        END DO
        END DO
    END DO
    

    DO k = 1, nwaves
        DO j = 1-nb, jhi+nb     
        DO i = 1-nb+1, ihi+nb-2
        
        
            IF(ABS(ck(i,j,k)) .GT. ZERO) THEN
                gammak = (gk(i+1,j,k) - gk(i,j,k)) / ck(i,j,k)
            ELSE
                gammak = 0.d0
            END IF
            
            chi = dtdx*eigenvalues(i,j,k) + gammak
       
            ! Qk_i+1/2 (K98 eqn. 2.31)
            IF(ABS(chi) .LT. (2.D0*eps(k))) THEN
                Qk = 0.25D0 * (chi**2 / eps(k)) + eps(k)
            ELSE
                Qk = ABS(chi)
            END IF   
            
            ! betak_i+1/2 (K98 eqn. 2.26)
            betak = Qk * ck(i,j,k) - (gk(i,j,k) + gk(i+1,j,k))  
            
            ! accumulate TVD flux contribution from each wave  (K98 eqn. 2.25)
            Flux(i,j,1,fcomp) = Flux(i,j,1,fcomp) - betak * Rk(i,j,1,k) 
            Flux(i,j,2,fcomp) = Flux(i,j,2,fcomp) - betak * Rk(i,j,2,k) 
            Flux(i,j,3,fcomp) = Flux(i,j,3,fcomp) - betak * Rk(i,j,3,k) 
            Flux(i,j,4,fcomp) = Flux(i,j,4,fcomp) - betak * Rk(i,j,4,k) 
            Flux(i,j,6,fcomp) = Flux(i,j,6,fcomp) - betak * Rk(i,j,5,k) 
            Flux(i,j,7,fcomp) = Flux(i,j,7,fcomp) - betak * Rk(i,j,6,k) 
            
        END DO    
        END DO
    END DO

    ! add remaining term in the Roe TVD flux
    DO j = 1-nb, jhi+nb     
    DO i = 1-nb+1, ihi+nb-2
        
        Flux(i,j,1,fcomp) = 0.5d0 * ( dxdt * Flux(i,j,1,fcomp) + F(i,j,1) + F(i+1,j,1) )
        Flux(i,j,2,fcomp) = 0.5d0 * ( dxdt * Flux(i,j,2,fcomp) + F(i,j,2) + F(i+1,j,2) )
        Flux(i,j,3,fcomp) = 0.5d0 * ( dxdt * Flux(i,j,3,fcomp) + F(i,j,3) + F(i+1,j,3) )
        Flux(i,j,4,fcomp) = 0.5d0 * ( dxdt * Flux(i,j,4,fcomp) + F(i,j,4) + F(i+1,j,4) )
        Flux(i,j,5,fcomp) = 0.d0
        Flux(i,j,6,fcomp) = 0.5d0 * ( dxdt * Flux(i,j,6,fcomp) + F(i,j,5) + F(i+1,j,5) )
        Flux(i,j,7,fcomp) = 0.5d0 * ( dxdt * Flux(i,j,7,fcomp) + F(i,j,6) + F(i+1,j,6) )   
        
    END DO
    END DO
    
    
END SUBROUTINE compute_TVD_fluxes_1d



! final (unsplit) finite volume update of the state vector (excluding Bx and By, which require CT update)
SUBROUTINE update_state_vector_2d()

    REAL*8 :: dtdx
    INTEGER :: i, j
    
    dtdx = dt / dx
    
    ! update the conserved variables (density, momentum and Bz only)
    DO j = 1-nb+2, ny+nb-2     
        DO i = 1-nb+2, nx+nb-2
        
            q_2d(i,j,1) = q_2d(i,j,1) - dtdx * ( (Fs(i,j,1,1) - Fs(i-1,j,1,1)) + (Fs(i,j,1,2) - Fs(i,j-1,1,2)) )
            q_2d(i,j,2) = q_2d(i,j,2) - dtdx * ( (Fs(i,j,2,1) - Fs(i-1,j,2,1)) + (Fs(i,j,2,2) - Fs(i,j-1,2,2)) )
            q_2d(i,j,3) = q_2d(i,j,3) - dtdx * ( (Fs(i,j,3,1) - Fs(i-1,j,3,1)) + (Fs(i,j,3,2) - Fs(i,j-1,3,2)) )
            q_2d(i,j,4) = q_2d(i,j,4) - dtdx * ( (Fs(i,j,4,1) - Fs(i-1,j,4,1)) + (Fs(i,j,4,2) - Fs(i,j-1,4,2)) )
            q_2d(i,j,5) = q_2d(i,j,5) - dtdx * ( (Fs(i,j,5,1) - Fs(i-1,j,5,1)) + (Fs(i,j,5,2) - Fs(i,j-1,5,2)) )   !! Get rid of these two if
            q_2d(i,j,6) = q_2d(i,j,6) - dtdx * ( (Fs(i,j,6,1) - Fs(i-1,j,6,1)) + (Fs(i,j,6,2) - Fs(i,j-1,6,2)) )   !! doing CT updates.
            q_2d(i,j,7) = q_2d(i,j,7) - dtdx * ( (Fs(i,j,7,1) - Fs(i-1,j,7,1)) + (Fs(i,j,7,2) - Fs(i,j-1,7,2)) )
            
        END DO
    END DO

    ! check for unphysical densities and apply protections accordingly
    DO j = -1, ny+1
        DO i = -1, nx+1
            IF(q_2d(i,j,1) .LT. ZERO) THEN
                PRINT*,'Negative density detected in Cell # ',i,j,', RHO= ',q_2d(i,j,1)
                STOP
            END IF  
        END DO    
    END DO


END SUBROUTINE update_state_vector_2d


! compute maximum speed (need this for calculating the time step size)
SUBROUTINE get_max_speed(max_speed)

    REAL*8, INTENT(INOUT) :: max_speed
    INTEGER :: i, j
    REAL*8 :: rho2, isrho, vx2, vy2, bx2, by2, bz2, bsqr, asqr, btsqr, cfx, cfy, tmp
    
    max_speed = 0.d0

    ! find the maximum fast mode speed on the grid (the fast mode speed at a cell-interface
    ! is computed using the roe-averaged state at that interface)
    DO j = 1-nb+2, ny+nb-3    
        DO i = 1-nb+2, nx+nb-3 ! no need to include boundary cells for this
            
            rho2 = 0.5d0 * ( q_2d(i,j,1) + q_2d(i+1,j,1) )
            isrho = 1.D0 / SQRT(rho2)
            vx2  = 0.5d0 * ( q_2d(i,j,2)/q_2d(i,j,1) + q_2d(i+1,j,2)/q_2d(i+1,j,1) )
            vy2  = 0.5d0 * ( q_2d(i,j,3)/q_2d(i,j,1) + q_2d(i+1,j,3)/q_2d(i+1,j,1) )
            bx2  = 0.5d0 * ( q_2d(i,j,5) + q_2d(i+1,j,5) ) * isrho
            by2  = 0.5d0 * ( q_2d(i,j,6) + q_2d(i+1,j,6) ) * isrho
            bz2  = 0.5d0 * ( q_2d(i,j,7) + q_2d(i+1,j,7) ) * isrho
        
            btsqr =  by2**2 + bz2**2 
            bsqr = bx2**2 + btsqr 
            asqr = sound_speed**2
            
            tmp = SQRT( (asqr - bsqr)**2 + 4.d0 * asqr * btsqr )
            cfx = SQRT( MAX(0.D0, 0.5d0 * (asqr + bsqr + tmp) ))  ! fast mode speed at x-interface


            rho2 = 0.5d0 * ( q_2d(i,j,1) + q_2d(i,j+1,1) )
            isrho = 1.D0 / SQRT(rho2)
            vx2  = 0.5d0 * ( q_2d(i,j,2)/q_2d(i,j,1) + q_2d(i,j+1,2)/q_2d(i,j+1,1) )
            vy2  = 0.5d0 * ( q_2d(i,j,3)/q_2d(i,j,1) + q_2d(i,j+1,3)/q_2d(i,j+1,1) )
            bx2  = 0.5d0 * ( q_2d(i,j,5) + q_2d(i,j+1,5) ) * isrho
            by2  = 0.5d0 * ( q_2d(i,j,6) + q_2d(i,j+1,6) ) * isrho
            bz2  = 0.5d0 * ( q_2d(i,j,7) + q_2d(i,j+1,7) ) * isrho
        
            btsqr =  by2**2 + bz2**2 
            bsqr  = bx2**2 + btsqr 
            asqr  = sound_speed**2
            
            tmp = SQRT( (asqr - bsqr)**2 + 4.d0 * asqr * btsqr )
            cfy = SQRT( MAX(0.D0, 0.5d0 * (asqr + bsqr + tmp) ))  ! fast mode speed at y-interface

            vx2 = q_2d(i,j,2)/q_2d(i,j,1)
            vy2 = q_2d(i,j,3)/q_2d(i,j,1)

            max_speed = MAX(max_speed, ABS(vx2) + cfx, ABS(vy2) + cfy)
            
        END DO
    END DO


    IF(max_speed .LT. ZERO) THEN
        PRINT*,'ERROR!! ZERO MAX SPEED!'
        STOP
    END IF


END SUBROUTINE get_max_speed



END MODULE Iso_TVD_Solver_mod