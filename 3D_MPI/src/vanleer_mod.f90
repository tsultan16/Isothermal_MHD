! 3D Isothermal MHD 2nd order Van Leer Solver (a.k.a. MUSCL-HANCOCK)
! References: Stone 'et al., New Astr., 14, 139-148, 2009  (S09)
!             Mignone, J. Comp., 225, 1427-1441, 2007 (M07)
!             Stone 'et al., ApJ. S., 178, 1, 2008  (S08)
!             Woodward+Colella, J.Comp, 54, 174-201, 1984 (WC)
!
!------------------------------------------------
! Outline of 3D Van Leer + CT Algorithm from S09:
!------------------------------------------------
! 
! The goal is to perform a directionally unsplit update of the Isothermal MHD state variables and
! maintain div(B) = 0 to machine precision.
!
! ****************
! * CT Algorithm * 
! ****************
! Div(B) = 0 is maintained by evolving area-averaged transverse magnetic fields 
! defined at each cell-faces via circulation of longitudinal electric fields 
! at the bounding cell edges, or cell-corners in 2D. (This is just the area integral
! version of Faraday's Law.). Details on how to compute the cell edge/corner electric 
! fields are given in (S18).
!
! 
! **********************
! * Van Leer Algorithm * 
! **********************
!
! The 3D update at each time step is performed simulatenously in both 
! x,y,z directions using a predictor-corrector approach.
!
! 1) Predictor Stage: We first compute first-order interface fluxes (and cell corner emf) using a Riemann solver.  
!                     The cell-center hydrodynamic variable are then updated by a half time-step (unsplit) using 
!                     these fluxes. The cell-face B-field undergoes a half time step CT update.
!
!
! 2) Corrector Stage: The MHD state variables and cell-face magnetic fields obtained from 
!                     the half-update in the predictor stage are then used to compute second
!                     order "interface states" reconstructed from a linear (or quadratic) interpolation of  
!                     cell-center primitive (or characteristic) variables. A slope limiter is applied during the 
!                     interpolation to ensure the TVD criterion. These interface states are then 
!                     used to compute the "corrector" cell-interface fluxes (and cell-corner electric fields)
!                     using the same Riemann solver. 
!                     These "corrector" fluxes and electric fields are then used to update 
!                     the cell-center hydrodynamic variables and cell-face B fields over the full time-step.
!
!
!
! ***************************************************************************************
! * Negative Density Protection Algorithm (only needed when using Roe's Riemann solver) *
! ***************************************************************************************
! When enabled, we concurrently compute "HLLE" fluxes. Each time the state variables get updated,
! we check for negative densities. If we detect a cell containing a negative density, we reupdate
! the state variables in that cell and in it's immediate neighbors using the HLLE fluxes at the 
! boundaries of the bad cell. Computing the HLLE fluxes adds an extra ~10% cost to overall solve
! per time-step. 
!-------------------------------------------------------------------------------------------------------------


MODULE VanLeer_Solver_mod

USE constants_mod
USE grid_data_mod
USE turb_driver_mod

IMPLICIT NONE


!*****************************************************************************************************************************!
!                                                     LOCAL VARIABLES                                                         !
!*****************************************************************************************************************************!


! Solver Parameters
INTEGER, PARAMETER :: nwaves = 6
REAL*8, PARAMETER  :: ZERO = 1.D-20           
REAL*8, PARAMETER :: MIN_DENS = 1.D-20        ! density floor
REAL*8, PARAMETER :: MC_beta = 1.01d0         ! MC limiter beta constant (1 <= beta <= 2), beta = 1 is equivalent to minmod and beta = 2 is MC
REAL*8 :: delta_h(6) = (/ 0.3d0, 0.d0, 0.3d0, 0.3d0, 0.d0, 0.3d0 /)  ! delta parameters for Harten's "entropy fix" for all six wave modes
REAL*8 :: steepen_var(6) = (/ 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0 /) ! parameters for specifying which primitive variable will be speened at a discontinuiuty, set this to either 1 or 0 for the desired variable

LOGICAL, PARAMETER :: flux_protection = .FALSE.          ! switch for turning on protection fluxes
LOGICAL, PARAMETER :: CT_off = .FALSE.                   ! this is can be used to switch off CT for testing purposes
LOGICAL, PARAMETER :: piecewise_parabolic = .TRUE.       ! switch for third-order spatial reconstruction for interface states
LOGICAL, PARAMETER :: primitive_reconstruction = .FALSE. ! switch for using primitive variable reconstruction rather than the default characteristic variable reconstruction 
                                                         ! (for piecewise linear, primitive variable reconstruction seems to be adequate. For piecewise parabolic, characteristic variable reconstruction is necessary)
LOGICAL, PARAMETER :: discontinuiuty_steepener = .FALSE. ! switch for discontinuiuty steepener

! Solver Local variables 
REAL*8, ALLOCATABLE :: qintL(:,:,:), qintR(:,:,:),   & ! state-vector 2D work arrays for interface-states
                       w(:,:,:), wb(:,:),            & ! primitive variables 2D work arrays
                       wL(:,:,:), wR(:,:,:),         & ! "
                       qwork_3d(:,:,:,:),        & ! state-vector 3D work array
                       q_3d_prot(:,:,:,:),       & ! storage for state vector protection backup
                       bface_work(:,:,:,:),      & ! cell-face magnetic field work array 
                       Fs(:,:,:,:,:),            & ! time-averaged flux vector
                       Fs_prot(:,:,:,:,:),       & ! storage for protection (HLLD) fluxes
                       Fs_2d(:,:,:),             & ! flux buffer for rotations
                       Fs_2d_prot(:,:,:),        & ! 2d array for protection flux storage
                       qavg(:,:,:),              & ! Roe-averaged state    
                       Rk(:,:,:,:),              & ! Roe-averaged right eigenvectors
                       ck(:,:,:),                & ! characteristic variables (i.e. projection of state vector onto the left eigenvectors)   
                       eigenvalues(:,:,:),       & ! Roe-averaged eigenvalues
                       wave_speeds(:,:,:),       & ! Roe-averaged wave speeds
                       emf_3d(:,:,:,:,:),        & ! cell-center reference and cell-face electric field z-component
                       emf_3d_prot(:,:,:,:,:),   & ! storage protection emf
                       emf_corner_prot(:,:,:,:), &  ! cell-corner CT electric field
                       Rw(:,:,:,:), Lw(:,:,:,:)  ! right and left-eigenvector matrices of the primitive variable eigensystem


LOGICAL, ALLOCATABLE :: protection_flags(:,:,:)


REAL*8, PARAMETER :: eps(6) = (/ 0.2, 0.0, 0.2, 0.2, 0.0 ,0.2/) ! wave dissipation constants

REAL*8 :: dt, dx, smax, delv_max

INTEGER :: nmax

!*****************************************************************************************************************************!
!                                                     SUBROUTINES                                                             !
!*****************************************************************************************************************************!

CONTAINS



SUBROUTINE init_solver_vanleer()
    
    ! make sure we have at least 5 boundary cells
    IF(nb .LT. 5) THEN
        PRINT*,'ERROR! Need nb >= 5 for Van Leer solver.'
        STOP
    END IF
    
    ! PPM requires an extra ghost cell layer
    IF(piecewise_parabolic .AND. nb .LT. 6) THEN
        PRINT*,'ERROR! Need nb >= 6 for third order spatial reconstruction of interfcace states.'
        STOP
    END IF

    ! make sure we are reconstructing in characteristic variables if using PPM
    IF(piecewise_parabolic .AND. primitive_reconstruction) THEN
        PRINT*,'WARNING! Use characteristic variable reconstruction for PPM!'
        !STOP
    END IF

    
    IF(myrank .EQ. 0) THEN
        PRINT*,''
        PRINT*,'Reconstruction in primitive variables: ', primitive_reconstruction
        PRINT*,'Using piecewise parabolic reconstruction: ', piecewise_parabolic
        PRINT*,''
    END IF
    
    CALL allocate_local_variables()


END SUBROUTINE init_solver_vanleer


SUBROUTINE destroy_solver_vanleer()

   CALL deallocate_local_variables()
   
END SUBROUTINE destroy_solver_vanleer


! allocate local variables/work arrays
SUBROUTINE allocate_local_variables()

    INTEGER :: total_grid_memory
    REAL*8 :: mbytes

    nmax = MAX(nx*nranks_x,ny*nranks_y,nz*nranks_z)

    ! set cell size  
    dx = 1.d0 / DBLE(nmax)
    
    ALLOCATE(qintL(1-nb:nmax+nb,1-nb:nmax+nb,7))
    ALLOCATE(qintR(1-nb:nmax+nb,1-nb:nmax+nb,7))
    ALLOCATE(w(1-nb:nmax+nb,1-nb:nmax+nb,6)) 
    ALLOCATE(wb(1-nb:nmax+nb,1-nb:nmax+nb)) 
    ALLOCATE(qwork_3d(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,7))
    ALLOCATE(bface_work(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb, 3))
    ALLOCATE(Fs(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,7,3))
    ALLOCATE(Fs_2d(1-nb:nmax+nb,1-nb:nmax+nb,7))
    ALLOCATE(emf_3d(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,3,3))

    IF(.NOT. primitive_reconstruction) THEN
        ALLOCATE(Rw(1-nb:nmax+nb,1-nb:nmax+nb,6,6))
        ALLOCATE(Lw(1-nb:nmax+nb,1-nb:nmax+nb,6,6))
    END IF   

    IF(piecewise_parabolic) ALLOCATE(wL(1-nb:nmax+nb,1-nb:nmax+nb,6), wR(1-nb:nmax+nb,1-nb:nmax+nb,6))


    IF(riemann_solver_type .EQ. 1) THEN
        ALLOCATE(qavg(1-nb:nmax+nb,1-nb:nmax+nb,7))
        ALLOCATE(Rk(1-nb:nmax+nb,1-nb:nmax+nb,nwaves, nwaves))
        ALLOCATE(ck(1-nb:nmax+nb,1-nb:nmax+nb,nwaves))
        ALLOCATE(eigenvalues(1-nb:nmax+nb,1-nb:nmax+nb,nwaves))
        ALLOCATE(wave_speeds(1-nb:nmax+nb,1-nb:nmax+nb,3))
    END IF

    IF(flux_protection) THEN

        ALLOCATE(q_3d_prot(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,7))
        ALLOCATE(Fs_prot(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,7,3))
        ALLOCATE(Fs_2d_prot(1-nb:nmax+nb,1-nb:nmax+nb,7))
        ALLOCATE(emf_3d_prot(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,3,3))
        ALLOCATE(emf_corner_prot(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,3))
        ALLOCATE(protection_flags(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb))

    END IF


    total_grid_memory = SIZEOF(qintL) + SIZEOF(qintR) +  SIZEOF(w) + SIZEOF(wb) +  SIZEOF(qwork_3d) + SIZEOF(bface_work) + &
                        SIZEOF(Fs) + SIZEOF(Fs_2d)  +  SIZEOF(emf_3d) 
    
    IF(.NOT. primitive_reconstruction) total_grid_memory = total_grid_memory + SIZEOF(Rw) + SIZEOF(Lw)
    
    IF(piecewise_parabolic) total_grid_memory = total_grid_memory + SIZEOF(wL) + SIZEOF(wR)
    
    IF(riemann_solver_type .EQ. 1) THEN
        total_grid_memory =  total_grid_memory + SIZEOF(qavg) + SIZEOF(Rk) + SIZEOF(ck) + SIZEOF(eigenvalues) + SIZEOF(wave_speeds)
    END IF

    IF(flux_protection) THEN
        total_grid_memory =  total_grid_memory + SIZEOF(q_3d_prot) + SIZEOF(Fs_prot) + &
                             SIZEOF(Fs_2d_prot) + SIZEOF(emf_3d_prot) + SIZEOF(emf_corner_prot) + SIZEOF(protection_flags)
    END IF

    mbytes = DBLE(total_grid_memory)*1.d-6    

    IF(myrank .EQ. 0) THEN
    
    PRINT*,''
    PRINT*,'#####################################################################'
    PRINT*,'Total Memory for MHD solver work array storage per MPI rank (Mb) = ', mbytes                   
    PRINT*,'#####################################################################'
    PRINT*,''

    END IF
    

END SUBROUTINE allocate_local_variables


SUBROUTINE deallocate_local_variables()


    DEALLOCATE(qintL, qintR, qwork_3d, bface_work)
    DEALLOCATE(w, wb) 
    DEALLOCATE(Fs)
    DEALLOCATE(emf_3d)
    DEALLOCATE(Fs_2d)

    IF(.NOT. primitive_reconstruction) DEALLOCATE(Rw, Lw)
    
    IF(piecewise_parabolic) DEALLOCATE(wL, wR)
    
    IF(riemann_solver_type .EQ. 1) THEN
        DEALLOCATE(qavg, Rk, ck, eigenvalues, wave_speeds)
    END IF

    IF(flux_protection) THEN 
        DEALLOCATE(q_3d_prot, Fs_prot, Fs_2d_prot)    
        DEALLOCATE(emf_3d_prot, emf_corner_prot)
        DEALLOCATE(protection_flags)
    END IF
    
    
END SUBROUTINE deallocate_local_variables


!*****************************************************************************************************************************!

! Top-level Isothermal TVD MHD solve routine (using CTU + CT scheme)
SUBROUTINE vanleer_solve_3d(dt_in, drive_next)

    REAL*8, INTENT(INOUT) :: dt_in
    LOGICAL, INTENT(IN) :: drive_next
    REAL*8 :: dtdx
    INTEGER :: i, j, offset  

    dt = dt_in
    
    ! First save initial state in work arrays
    qwork_3d = q_3d
    bface_work = bface_3d

    ! clear flux arrays (may not be necessary)
    Fs = 0.d0
    Fs_prot = 0.d0

    ! clear electric field array (may not be necessary)
    emf_3d = 0.d0
    emf_corner = 0.d0  
 
    ! reset max wave speed
    smax = 0.d0
    delv_max = 0.d0
    
    ! send initial density to turbulence driver
    IF(drive_turbulence .AND. drive_next) CALL get_rho_turb()
    
    offset = 0
    
    !*****************************************************************************************************************
    ! Step 1: Compute first order predictor cell-interafce fluxes along x, y and z directions using 1D Riemann solver
    !*****************************************************************************************************************

    IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing predictor TVD x-fluxes.'
    
    CALL compute_xfluxes_3d(1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset, .FALSE.) 

    IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing predictor TVD y-fluxes.'
    
    CALL compute_yfluxes_3d(1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset, .FALSE.) 

    IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing predictor TVD z-fluxes.'
    
    CALL compute_zfluxes_3d(1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset, .FALSE.) 
            
    !***********************************************************************************
    ! Step 2: Compute predictor cell-centered reference and cell-corner electric fields  
    !***********************************************************************************    
    
    IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing predictor CT cell center reference electric field.'

    CALL compute_cell_center_emf_3d(1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset)  

    IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing predictor CT cell corner electric field.'

    CALL compute_cell_corner_emf_3d(emf_3d, emf_corner, 1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset)
    IF(flux_protection) CALL compute_cell_corner_emf_3d(emf_3d_prot, emf_corner_prot, 1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset)

    !******************************************************************
    ! Step 3: Apply 1/2 dt update to cell-center hydrodynamic variables
    !******************************************************************

    dtdx = 0.5d0 * dt / dx
   
    IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing hydro state half update.'
    
    CALL update_hydro_state_3d(dtdx, 1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset)

    !*************************************************************************************
    ! Step 3 & 4: Apply 1/2 dt CT update to cell-face magnetic fields and interpolate
    !             to get cell-center magnetic fields
    !************************************************************************************
    
    IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing CT half update.'

    IF(.NOT. CT_off) CALL update_CT_3d(dtdx, 1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset) 
   
    offset = 2
      
    !************************************************************************************************
    ! Step 5 & 6: Apply piecewise linear reconstruction to compute interface-states at half time-step
    !             then use these to compute the corrector interface fluxes
    !********************************************************************************************
       
    ! compute corrector x-fluxes    
    IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing corrector fluxes in y direction.'
    
    CALL compute_xfluxes_3d(1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset, .TRUE.)  
         
    ! compute corrector y-fluxes     
    IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing corrector fluxes in y direction.'

    CALL compute_yfluxes_3d(1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset, .TRUE.)      

    ! compute corrector z-fluxes
    IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing corrector fluxes in z direction.'

    CALL compute_zfluxes_3d(1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset, .TRUE.) 

    IF(.NOT. piecewise_parabolic) THEN
        offset = 3 
    ELSE
        offset = 4
    END IF    

    !***************************************************************************************
    ! Step 6 & 7: Compute corrector cell-centered reference and cell-corner electric fields  
    !***************************************************************************************
    
    IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing corrector CT cell center reference electric field.'

    CALL compute_cell_center_emf_3d(1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset)

    IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing corrector CT cell corner electric field.'

    CALL compute_cell_corner_emf_3d(emf_3d, emf_corner, 1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset)
    IF(flux_protection) CALL compute_cell_corner_emf_3d(emf_3d_prot, emf_corner_prot, 1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset)

    !***********************************************************************************
    ! Step 9: Apply final full-dt update to cell-center density, momentum and cell-face  
    !         magnetic fields using the corrector fluxes and cell-corner electic fields
    !*********************************************************************************** 

    ! Restore intial state 
    q_3d = qwork_3d
    bface_3d = bface_work
        
    dtdx = dt / dx

    IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing final unsplit update of state variables and CT.'

    CALL update_hydro_state_3d(dtdx, 1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset)

    IF(.NOT. CT_off) CALL update_CT_3d(dtdx, 1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset) 
    
    !***********************************************************************
    ! Step 10: Add driving force terms for turbulence 
    !***********************************************************************
    
    IF(drive_turbulence .AND. drive_next) THEN
        CALL add_random_force(dt, delv_re, delv_max)
    END IF   
    
    IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Update for time step completed.'

    
END SUBROUTINE vanleer_solve_3d


SUBROUTINE compute_xfluxes_3d(ilow, ihi, jlow, jhi, klow, khi, corrector_stage) 

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi, klow, khi
    LOGICAL, INTENT(IN) :: corrector_stage
    INTEGER :: i, j, k, offset
    REAL*8  :: irho 
 
    ! Loop over z planes    
    DO k = klow-nb, khi+nb
 
        IF(print_debug .AND. myrank .EQ. 0) PRINT*,' Z plane index,  k = ', k
        
        ! Compute interface states. (For predictor stage, the interface states needed for first order fluxes
        ! are just cell-center values) 
        IF(.NOT. corrector_stage) THEN
            
            DO j = jlow-nb, jhi+nb 
                DO i = ilow-nb, ihi+nb-1        

                    ! left state
                    qintL(i,j,1) = q_3d(i,j,k,1)
                    qintL(i,j,2) = q_3d(i,j,k,2)
                    qintL(i,j,3) = q_3d(i,j,k,3)
                    qintL(i,j,4) = q_3d(i,j,k,4)
                    qintL(i,j,5) = bface_3d(i,j,k,1) 
                    qintL(i,j,6) = q_3d(i,j,k,6)
                    qintL(i,j,7) = q_3d(i,j,k,7)

                    ! right state
                    qintR(i,j,1) = q_3d(i+1,j,k,1)
                    qintR(i,j,2) = q_3d(i+1,j,k,2)
                    qintR(i,j,3) = q_3d(i+1,j,k,3)
                    qintR(i,j,4) = q_3d(i+1,j,k,4)
                    qintR(i,j,5) = bface_3d(i,j,k,1) 
                    qintR(i,j,6) = q_3d(i+1,j,k,6)
                    qintR(i,j,7) = q_3d(i+1,j,k,7)
                    
                END DO
            END DO
            
            offset = 0 

        ELSE
        
            ! first compute cell-center primitive variables
            DO j = jlow-nb, jhi+nb 
                DO i = ilow-nb, ihi+nb                     
            
                    irho = 1.d0 / q_3d(i,j,k,1)
          
                    w(i,j,1) = q_3d(i,j,k,1)
                    w(i,j,2) = q_3d(i,j,k,2) * irho
                    w(i,j,3) = q_3d(i,j,k,3) * irho
                    w(i,j,4) = q_3d(i,j,k,4) * irho
                    wb(i,j) = bface_3d(i,j,k,1)
                    w(i,j,5) = q_3d(i,j,k,6) 
                    w(i,j,6) = q_3d(i,j,k,7) 
                    
                    ! Note: The normal component of the magnetic field is already defined at the cell interface
                    ! so will not require interpolation
            
                END DO
            END DO
                    
            ! Now interpolate to get cell-iterface primitive variables and reconstruct cell-interface states from that 
            ! (for corrector stage, generating higher-order fluxes requires that we interpolate from cell-center to the interface values)
            IF(.NOT. primitive_reconstruction) THEN
                CALL compute_primitive_eigenvectors(ilow, ihi, jlow, jhi)
                CALL compute_interface_states_char(ilow, ihi, jlow, jhi)
            ELSE
                CALL compute_interface_states_prim(ilow, ihi, jlow, jhi)
            END IF        
            
            IF(.NOT. piecewise_parabolic) THEN
                offset = 1 
            ELSE
                offset = 2
            END IF       
                
        END IF

       
        IF(riemann_solver_type .EQ. 3) THEN
        
            CALL compute_HLLE_fluxes_1d(Fs_2d, ilow+offset, ihi-offset, jlow, jhi)
        
        ELSE IF(riemann_solver_type .EQ. 2) THEN
        
            CALL compute_HLLD_fluxes_1d(Fs_2d, ilow+offset, ihi-offset, jlow, jhi)
        
        
        ELSE IF(riemann_solver_type .EQ. 1) THEN
 
 
            IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing Roe average state.' 
    
            ! now sweep through strips along x-direction and compute x-interface fluxes   
            CALL compute_roe_avg_state_1d(ilow+offset, ihi-offset, jlow, jhi)
    
            IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing eigenvalues.' 
    
            CALL compute_roe_eigenvalues_1d(ilow+offset, ihi-offset, jlow, jhi)
            
            IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing eigenvectors.'
        
            CALL compute_roe_eigenvectors_1d(ilow+offset, ihi-offset, jlow, jhi)

            IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing x-fluxes.'

            CALL compute_ROE_fluxes_1d(Fs_2d, ilow+offset, ihi-offset, jlow, jhi)

            IF(flux_protection) CALL compute_HLLE_fluxes_1d(Fs_2d_prot, ilow+offset, ihi-offset, jlow, jhi)

        END IF
                
        
        DO j = jlow-nb, jhi+nb 
            DO i = ilow+offset-nb, ihi-offset+nb-1   
                Fs(i,j,k,:,1) = Fs_2d(i,j,:)
            END DO
        END DO
            
        ! this is a good place to store the cell-face transverse electric fields
        DO j = jlow-nb, jhi+nb
            DO i = ilow+offset-nb, ihi-offset+nb-1   
                emf_3d(i,j,k,3,1) = -Fs(i,j,k,6,1)   ! Flux-x(By)_i+1/2,j,k = -Ez_i+1/2,j,k
                emf_3d(i,j,k,2,1) = Fs(i,j,k,7,1)    ! Flux-x(Bz)_i+1/2,j,k = Ey_i+1/2,j,k
            END DO
        END DO    
    
    
        IF(flux_protection) THEN
    
            DO j = jlow-nb, jhi+nb 
                DO i = ilow-nb, ihi+nb 
                    Fs_prot(i,j,k,:,1) = Fs_2d_prot(i,j,:)
                END DO
            END DO
        
            ! this is a good place to store the cell-face transverse electric fields
            DO j = jlow-nb, jhi+nb
                DO i = ilow+offset-nb, ihi-offset+nb-1     
                    emf_3d_prot(i,j,k,3,1) = -Fs_prot(i,j,k,6,1)   ! Flux-x(By)_i+1/2,j,k = -Ez_i+1/2,j,k 
                    emf_3d_prot(i,j,k,2,1) = Fs_prot(i,j,k,7,1)    ! Flux-x(Bz)_i+1/2,j,k = Ey_i+1/2,j,k        
                END DO
            END DO    
    
        END IF
    
    
    END DO
    

END SUBROUTINE compute_xfluxes_3d 


SUBROUTINE compute_yfluxes_3d(ilow, ihi, jlow, jhi, klow, khi, corrector_stage) 

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi, klow, khi
    LOGICAL, INTENT(IN) :: corrector_stage
    INTEGER :: i, j, k, offset
    REAL*8 :: irho

    ! loop over z planes
    DO k = klow-nb, khi+nb
    
        IF(print_debug .AND. myrank .EQ. 0) PRINT*,' Z plane index,  k = ', k
   

        ! Compute interface states. (For predictor stage, the interface states needed for first order fluxes
        ! are just cell-center values)      
        IF(.NOT. corrector_stage) THEN
            
            
            ! rotate the stave vector array so that fastest index runs along y-direction
            ! (also need to permute x->y->z components for momentum and B field)
            
            DO i = ilow-nb, ihi+nb 
                DO j = jlow-nb, jhi+nb-1
            
                    ! left state
                    qintL(j,i,1) = q_3d(i,j,k,1)
                    qintL(j,i,2) = q_3d(i,j,k,3)
                    qintL(j,i,3) = q_3d(i,j,k,4)
                    qintL(j,i,4) = q_3d(i,j,k,2)
                    qintL(j,i,5) = bface_3d(i,j,k,2) 
                    qintL(j,i,6) = q_3d(i,j,k,7)
                    qintL(j,i,7) = q_3d(i,j,k,5)
            
                    ! right state
                    qintR(j,i,1) = q_3d(i,j+1,k,1)
                    qintR(j,i,2) = q_3d(i,j+1,k,3)
                    qintR(j,i,3) = q_3d(i,j+1,k,4)
                    qintR(j,i,4) = q_3d(i,j+1,k,2)
                    qintR(j,i,5) = bface_3d(i,j,k,2) 
                    qintR(j,i,6) = q_3d(i,j+1,k,7)
                    qintR(j,i,7) = q_3d(i,j+1,k,5)
            
            
                END DO
            END DO
     
            offset = 0
     
        ELSE
        
            ! first compute cell-center primitive variables,
            ! also rotate so that fastest index runs along y-direction
            ! (also need to permute  x->y->z components for momentum and B field)
            DO i = ilow-nb, ihi+nb 
                DO j = jlow-nb, jhi+nb                     
            
                    irho = 1.d0 / q_3d(i,j,k,1)
                              
                    w(j,i,1) = q_3d(i,j,k,1)
                    w(j,i,2) = q_3d(i,j,k,3) * irho
                    w(j,i,3) = q_3d(i,j,k,4) * irho
                    w(j,i,4) = q_3d(i,j,k,2) * irho
                    wb(j,i) = bface_3d(i,j,k,2)
                    w(j,i,5) = q_3d(i,j,k,7)
                    w(j,i,6) = q_3d(i,j,k,5)
                    
                    
                    ! Note: The normal component of the magnetic field is already defined at the cell interface
                    ! so will not require interpolation
                               
                END DO
            END DO
            
            ! Now interpolate to get cell-iterface primitive variables and reconstruct cell-interface states from that 
            ! (for corrector stage, generating second order fluxes requires that we linearly interpolate from cell-center to the interface values)
            IF(.NOT. primitive_reconstruction) THEN
                CALL compute_primitive_eigenvectors(jlow, jhi, ilow, ihi)
                CALL compute_interface_states_char(jlow, jhi, ilow, ihi)
            ELSE
                CALL compute_interface_states_prim(jlow, jhi, ilow, ihi)
            END IF
                
            IF(.NOT. piecewise_parabolic) THEN
                offset = 1 
            ELSE
                offset = 2
            END IF       
                
        END IF
     
        
        IF(riemann_solver_type .EQ. 3) THEN
        
            CALL compute_HLLE_fluxes_1d(Fs_2d, jlow+offset, jhi-offset, ilow, ihi)
        
        ELSE IF(riemann_solver_type .EQ. 2) THEN
        
            CALL compute_HLLD_fluxes_1d(Fs_2d, jlow+offset, jhi-offset, ilow, ihi)
        
        ELSE IF(riemann_solver_type .EQ. 1) THEN
        
            IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing Roe average state.' 

            ! now sweep through strips along y-direction and compute y-interface fluxes   
            CALL compute_roe_avg_state_1d(jlow+offset, jhi-offset, ilow, ihi)
    
            IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing eigenvalues.'
     
            CALL compute_roe_eigenvalues_1d(jlow+offset, jhi-offset, ilow, ihi)
       
            IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing eigenvectors.'
    
            CALL compute_roe_eigenvectors_1d(jlow+offset, jhi-offset, ilow, ihi)

            IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing y-fluxes.'

            CALL compute_ROE_fluxes_1d(Fs_2d, jlow+offset, jhi-offset, ilow, ihi)

            IF(flux_protection) CALL compute_HLLE_fluxes_1d(Fs_2d_prot, jlow+offset, jhi-offset, ilow, ihi)

        END IF

       
        ! Rotate the y-fluxes in flux buffer and store them in main flux array 
        DO j = jlow+offset-nb, jhi-offset+nb-1
            DO i = ilow-nb, ihi+nb
        
                Fs(i,j,k,1,2) = Fs_2d(j,i,1)
                Fs(i,j,k,3,2) = Fs_2d(j,i,2)
                Fs(i,j,k,4,2) = Fs_2d(j,i,3)
                Fs(i,j,k,2,2) = Fs_2d(j,i,4)
                Fs(i,j,k,6,2) = Fs_2d(j,i,5)
                Fs(i,j,k,7,2) = Fs_2d(j,i,6)         
                Fs(i,j,k,5,2) = Fs_2d(j,i,7)           
            
            END DO
        END DO
    
        ! this is a good place to store the cell-face transverse electric fields
        DO j = jlow+offset-nb, jhi-offset+nb-1
            DO i = ilow-nb, ihi+nb  
                emf_3d(i,j,k,1,1) = -Fs(i,j,k,7,2) ! Flux-y(Bz)_i,j+1/2,k = -Ex_i,j+1/2,k
                emf_3d(i,j,k,3,2) = Fs(i,j,k,5,2)  ! Flux-y(Bx)_i,j+1/2,k = Ez_i,j+1/2,k
            END DO
        END DO      
            
        IF(flux_protection) THEN
 
            ! Rotate the y-fluxes in flux buffer and store them in main flux array 
            DO j = jlow+offset-nb, jhi-offset+nb-1
                DO i = ilow-nb, ihi+nb
        
                    Fs_prot(i,j,k,1,2) = Fs_2d_prot(j,i,1)
                    Fs_prot(i,j,k,3,2) = Fs_2d_prot(j,i,2)
                    Fs_prot(i,j,k,4,2) = Fs_2d_prot(j,i,3)
                    Fs_prot(i,j,k,2,2) = Fs_2d_prot(j,i,4)
                    Fs_prot(i,j,k,6,2) = Fs_2d_prot(j,i,5)
                    Fs_prot(i,j,k,7,2) = Fs_2d_prot(j,i,6)         
                    Fs_prot(i,j,k,5,2) = Fs_2d_prot(j,i,7)           
                
                END DO
            END DO
    
            ! this is a good place to store the cell-face transverse electric fields
            DO j = jlow+offset-nb, jhi-offset+nb-1
                DO i = ilow-nb, ihi+nb  
                    emf_3d_prot(i,j,k,1,1) = -Fs_prot(i,j,k,7,2)   ! Flux-y(Bz)_i,j+1/2,k = -Ex_i,j+1/2,k
                    emf_3d_prot(i,j,k,3,2) = Fs_prot(i,j,k,5,2)  ! Flux-y(Bx)_i,j+1/2,k = -Ez_i,j+1/2,k
                END DO
            END DO                      
 
        END IF
    

    END DO

END SUBROUTINE compute_yfluxes_3d 


SUBROUTINE compute_zfluxes_3d(ilow, ihi, jlow, jhi, klow, khi, corrector_stage) 

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi, klow, khi
    LOGICAL, INTENT(IN) :: corrector_stage
    INTEGER :: i, j, k, offset
    REAL*8  :: irho 


    ! loop over y planes
    DO j = jlow-nb, jhi+nb  
    
        IF(print_debug .AND. myrank .EQ. 0) PRINT*,' Y plane index:  j = ', j
    
        ! Compute interface states. (For predictor stage, the interface states needed for first order fluxes
        ! are just cell-center values)      
        IF(.NOT. corrector_stage) THEN
            
            
            ! first rotate the stave vector array so that fastest index runs along z-direction 
            ! (also need to permute x->y->z components for momentum and B field)
            DO i = ilow-nb, ihi+nb 
                DO k = klow-nb, khi+nb-1
            
                    ! left state
                    qintL(k,i,1) = q_3d(i,j,k,1)
                    qintL(k,i,2) = q_3d(i,j,k,4)
                    qintL(k,i,3) = q_3d(i,j,k,2)
                    qintL(k,i,4) = q_3d(i,j,k,3)
                    qintL(k,i,5) = bface_3d(i,j,k,3) 
                    qintL(k,i,6) = q_3d(i,j,k,5)
                    qintL(k,i,7) = q_3d(i,j,k,6)
            
                    ! right state
                    qintR(k,i,1) = q_3d(i,j,k+1,1)
                    qintR(k,i,2) = q_3d(i,j,k+1,4)
                    qintR(k,i,3) = q_3d(i,j,k+1,2)
                    qintR(k,i,4) = q_3d(i,j,k+1,3)
                    qintR(k,i,5) = bface_3d(i,j,k,3)
                    qintR(k,i,6) = q_3d(i,j,k+1,5)
                    qintR(k,i,7) = q_3d(i,j,k+1,6)
            
                END DO
            END DO
            
            offset = 0
            
        ELSE
        
            ! first compute cell-center primitive variables,
            ! also rotate so that fastest index runs along z-direction
            ! (also need to permute  x->y->z components for momentum and B field)
            DO i = ilow-nb, ihi+nb 
                DO k = klow-nb, khi+nb                   
            
                    irho = 1.d0 / q_3d(i,j,k,1)
                              
                    w(k,i,1) = q_3d(i,j,k,1)
                    w(k,i,2) = q_3d(i,j,k,4) * irho
                    w(k,i,3) = q_3d(i,j,k,2) * irho
                    w(k,i,4) = q_3d(i,j,k,3) * irho
                    wb(k,i) = bface_3d(i,j,k,3) 
                    w(k,i,5) = q_3d(i,j,k,5)
                    w(k,i,6) = q_3d(i,j,k,6)
                    
                    ! Note: The normal component of the magnetic field is already defined at the cell interface
                    ! so will not require interpolation
                               
                END DO
            END DO
            
            ! Now interpolate to get cell-iterface primitive variables and reconstruct cell-interface states from that 
            ! (for corrector stage, generating second order fluxes requires that we linearly interpolate from cell-center to the interface values)
            
            IF(.NOT. primitive_reconstruction) THEN
                CALL compute_primitive_eigenvectors(klow, khi, ilow, ihi)
                CALL compute_interface_states_char(klow, khi, ilow, ihi)
            ELSE
                CALL compute_interface_states_prim(klow, khi, ilow, ihi)
            END IF            
                
            IF(.NOT. piecewise_parabolic) THEN
                offset = 1 
            ELSE
                offset = 2
            END IF       
        
        END IF
       
        
        IF(riemann_solver_type .EQ. 3) THEN
        
            CALL compute_HLLE_fluxes_1d(Fs_2d, klow+offset, khi-offset, ilow, ihi)
        
        ELSE IF(riemann_solver_type .EQ. 2) THEN
        
            CALL compute_HLLD_fluxes_1d(Fs_2d, klow+offset, khi-offset, ilow, ihi)
        
        ELSE IF(riemann_solver_type .EQ. 1) THEN
        
            IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing Roe average state.' 

            ! now sweep through strips along y-direction and compute y-interface fluxes   
            CALL compute_roe_avg_state_1d(klow+offset, khi-offset, ilow, ihi)
    
            IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing eigenvalues.'
     
            CALL compute_roe_eigenvalues_1d(klow+offset, khi-offset, ilow, ihi)
       
            IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing eigenvectors.'
    
            CALL compute_roe_eigenvectors_1d(klow+offset, khi-offset, ilow, ihi)

            IF(print_debug .AND. myrank .EQ. 0) PRINT*,'Computing z-fluxes.'

            CALL compute_ROE_fluxes_1d(Fs_2d, klow+offset, khi-offset, ilow, ihi)

            IF(flux_protection) CALL compute_HLLE_fluxes_1d(Fs_2d_prot, klow+offset, khi-offset, ilow, ihi)

        END IF

        
        ! Rotate the z-fluxes in flux buffer and store them in main flux array 
        DO k = klow+offset-nb, khi-offset+nb-1
            DO i = ilow-nb, ihi+nb
      
                Fs(i,j,k,1,3) = Fs_2d(k,i,1)
                Fs(i,j,k,4,3) = Fs_2d(k,i,2)
                Fs(i,j,k,2,3) = Fs_2d(k,i,3)
                Fs(i,j,k,3,3) = Fs_2d(k,i,4)
                Fs(i,j,k,7,3) = Fs_2d(k,i,5)
                Fs(i,j,k,5,3) = Fs_2d(k,i,6)         
                Fs(i,j,k,6,3) = Fs_2d(k,i,7)         
            
            END DO
        END DO
    
        ! this is a good place to store the cell-face transverse electric fields
        DO k = klow+offset-nb, khi-offset+nb-1
            DO i = ilow-nb, ihi+nb   
                emf_3d(i,j,k,2,2) = -Fs(i,j,k,5,3)   ! Flux-z(Bx)_i,j,k+1/2 = -Ey_i,j,k+1/2
                emf_3d(i,j,k,1,2) = Fs(i,j,k,6,3)    ! Flux-z(By)_i,j,k+1/2 = Ex_i,j,k+1/2              
            END DO
        END DO      
            
        IF(flux_protection) THEN

            ! Rotate the z-fluxes in flux buffer and store them in main flux array 
            DO k = klow+offset-nb, khi-offset+nb-1
                DO i = ilow-nb, ihi+nb
        
                    Fs_prot(i,j,k,1,3) = Fs_2d_prot(k,i,1)
                    Fs_prot(i,j,k,4,3) = Fs_2d_prot(k,i,2)
                    Fs_prot(i,j,k,2,3) = Fs_2d_prot(k,i,3)
                    Fs_prot(i,j,k,3,3) = Fs_2d_prot(k,i,4)
                    Fs_prot(i,j,k,7,3) = Fs_2d_prot(k,i,5)
                    Fs_prot(i,j,k,5,3) = Fs_2d_prot(k,i,6)         
                    Fs_prot(i,j,k,6,3) = Fs_2d_prot(k,i,7)         
            
                END DO
            END DO
    
            ! this is a good place to store the cell-face transverse electric fields
            DO k = klow+offset-nb, khi-offset+nb-1
                DO i = ilow-nb, ihi+nb   
                    emf_3d_prot(i,j,k,2,2) = -Fs_prot(i,j,k,5,3)   ! Flux-z(Bx)_i,j,k+1/2 = Ey_i,j,k+1/2
                    emf_3d_prot(i,j,k,1,2) = Fs_prot(i,j,k,6,3)    ! Flux-z(By)_i,j,k+1/2 = Ex_i,j,k+1/2              
                END DO
            END DO    

        END IF
    

    END DO

END SUBROUTINE compute_zfluxes_3d 


SUBROUTINE  compute_primitive_eigenvectors(ilow, ihi, jlow, jhi)

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi
    INTEGER :: i, j
    REAL*8 :: rho, irho, cs, ca, cf, bx, by, bz, asqr, bsqr, btsqr, rhos, irhos, tmp1, tmp2
    REAL*8 :: alphaf, alphas, betay, betaz, Cff, Css, Qf, Qs, Af, As, S, Nf, Ns 

   
    asqr = sound_speed**2

    ! matrix of right eigenvectors , ith column is the ith eigenvector  (see S08 eqn. A12)
    DO j = jlow-nb, jhi+nb 
        DO i = ilow-nb, ihi+nb

            ! compute the wave speeds
            rho = w(i,j,1)
            rhos = DSQRT(rho)
            irhos = 1.d0 / rhos
            bx = irhos * wb(i,j) 
            by = irhos * w(i,j,5) 
            bz = irhos * w(i,j,6) 
            bsqr  = bx**2 + by**2 + bz**2
            btsqr = by**2 + bz**2
            tmp1  = asqr + bsqr        
            tmp2  = DSQRT( (asqr - bsqr)**2 + 4.d0 * asqr * btsqr )

            cs = DSQRT( MAX(0.D0, 0.5d0 * (tmp1 - tmp2) ))  ! slow mode speed
            ca = ABS(bx)                                    ! alfven speed
            cf = DSQRT(0.5d0 * (tmp1 + tmp2))               ! fast mode speed

            ! compute auxiliary constants
            !IF(ABS(ca**2-sound_speed**2) .LT. ZERO .OR. ABS(cf**2-cs**2) .LT. ZERO) THEN
            IF(ABS(cf**2-cs**2) .LT. ZERO) THEN
                alphaf = 1.d0 
                alphas = 0.d0 
            ELSE
                alphaf = DSQRT(MAX((asqr - cs**2) / (cf**2 - cs**2), 0.d0))
                alphas = DSQRT(MAX((cf**2 - asqr) / (cf**2 - cs**2), 0.d0))
            END IF
            
            IF(btsqr .LT. ZERO) THEN
                betay = 1.D0 / DSQRT(2.D0)
                betaz = 1.D0 / DSQRT(2.D0)
            ELSE
                betay = by / DSQRT(btsqr)
                betaz = bz / DSQRT(btsqr)
            END IF
            
            S   = SIGN(1.d0, bx)
            Cff = cf * alphaf     
            Css = cs * alphas    
            Qf  = cf * alphaf * S 
            Qs  = cs * alphas * S
            Af  = sound_speed * alphaf * rhos
            As  = sound_speed * alphas * rhos

            ! compute right eigenvector columns
            Rw(i,j,1,1) = rho * alphaf
            Rw(i,j,2,1) = -Cff
            Rw(i,j,3,1) = Qs * betay
            Rw(i,j,4,1) = Qs * betaz
            Rw(i,j,5,1) = As * betay
            Rw(i,j,6,1) = As * betaz

            Rw(i,j,1,2) = 0.d0
            Rw(i,j,2,2) = 0.d0
            Rw(i,j,3,2) = -betaz
            Rw(i,j,4,2) = betay
            Rw(i,j,5,2) = -betaz * S * rhos
            Rw(i,j,6,2) = betay * S * rhos
            
            Rw(i,j,1,3) = rho * alphas
            Rw(i,j,2,3) = -Css
            Rw(i,j,3,3) = -Qf * betay
            Rw(i,j,4,3) = -Qf * betaz
            Rw(i,j,5,3) = -Af * betay
            Rw(i,j,6,3) = -Af * betaz
            
            Rw(i,j,1,4) = rho * alphas
            Rw(i,j,2,4) = Css
            Rw(i,j,3,4) = Qf * betay
            Rw(i,j,4,4) = Qf * betaz
            Rw(i,j,5,4) = -Af * betay
            Rw(i,j,6,4) = -Af * betaz
            
            Rw(i,j,1,5) = 0.d0
            Rw(i,j,2,5) = 0.d0
            Rw(i,j,3,5) = betaz
            Rw(i,j,4,5) = -betay
            Rw(i,j,5,5) = -betaz * S * rhos
            Rw(i,j,6,5) = betay * S * rhos
            
            Rw(i,j,1,6) = rho * alphaf
            Rw(i,j,2,6) = Cff
            Rw(i,j,3,6) = -Qs * betay
            Rw(i,j,4,6) = -Qs * betaz
            Rw(i,j,5,6) = As * betay
            Rw(i,j,6,6) = As * betaz
            
        END DO
    END DO
    
    
    ! matrix of left eigenvectors , ith row is the ith eigenvector  (see S08 eqn. A22)
    DO j = jlow-nb, jhi+nb 
        DO i = ilow-nb, ihi+nb

            ! compute the wave speeds
            rho = w(i,j,1)
            irho = 1.d0 / rho
            rhos = DSQRT(rho)
            irhos = 1.d0 / rhos
            bx = irhos * wb(i,j) 
            by = irhos * w(i,j,5) 
            bz = irhos * w(i,j,6) 
            bsqr  = bx**2 + by**2 + bz**2
            btsqr = by**2 + bz**2
            tmp1  = asqr + bsqr        
            tmp2  = DSQRT( (asqr - bsqr)**2 + 4.d0 * asqr * btsqr )

            cs = DSQRT( MAX(0.D0, 0.5d0 * (tmp1 - tmp2) ))  ! slow mode speed
            ca = ABS(bx)                                    ! alfven speed
            cf = DSQRT(0.5d0 * (tmp1 + tmp2))               ! fast mode speed

            ! compute auxiliary constants
            !IF(ABS(ca**2-sound_speed**2) .LT. ZERO .OR. ABS(cf**2-cs**2) .LT. ZERO) THEN
            IF(ABS(cf**2-cs**2) .LT. ZERO) THEN
                alphaf = 1.d0 
                alphas = 0.d0 
            ELSE
                alphaf = DSQRT(MAX((asqr - cs**2) / (cf**2 - cs**2), 0.d0))
                alphas = DSQRT(MAX((cf**2 - asqr) / (cf**2 - cs**2), 0.d0))
            END IF
            
            IF(btsqr .LT. ZERO) THEN
                betay = 1.D0 / DSQRT(2.D0)
                betaz = 1.D0 / DSQRT(2.D0)
            ELSE
                betay = by / DSQRT(btsqr)
                betaz = bz / DSQRT(btsqr)
            END IF
            
            S   = SIGN(1.d0, bx)
            Cff = cf * alphaf     
            Css = cs * alphas    
            Qf  = cf * alphaf * S 
            Qs  = cs * alphas * S
            Af  = sound_speed * alphaf * rhos
            As  = sound_speed * alphas * rhos
            Nf  = 0.5d0 / asqr  
            Ns  = Nf
            
            ! compute left eigenvector rows
            Lw(i,j,1,1) = Nf * alphaf * asqr * irho
            Lw(i,j,1,2) = -Nf * Cff
            Lw(i,j,1,3) = Nf * Qs * betay
            Lw(i,j,1,4) = Nf * Qs * betaz
            Lw(i,j,1,5) = Nf * As * betay * irho
            Lw(i,j,1,6) = Nf * As * betaz * irho
            
            Lw(i,j,2,1) = 0.d0
            Lw(i,j,2,2) = 0.d0
            Lw(i,j,2,3) = -0.5d0 * betaz
            Lw(i,j,2,4) = 0.5d0 * betay
            Lw(i,j,2,5) = -0.5d0 * betaz * S * irhos
            Lw(i,j,2,6) = 0.5d0 * betay * S * irhos
            
            Lw(i,j,3,1) = Ns * alphas * asqr * irho
            Lw(i,j,3,2) = -Ns * Css
            Lw(i,j,3,3) = -Ns * Qf * betay
            Lw(i,j,3,4) = -Ns * Qf * betaz
            Lw(i,j,3,5) = -Ns * Af * betay * irho
            Lw(i,j,3,6) = -Ns * Af * betaz * irho
            
            Lw(i,j,4,1) = Ns * alphas * asqr * irho
            Lw(i,j,4,2) = Ns * Css
            Lw(i,j,4,3) = Ns * Qf * betay
            Lw(i,j,4,4) = Ns * Qf * betaz
            Lw(i,j,4,5) = -Ns * Af * betay * irho
            Lw(i,j,4,6) = -Ns * Af * betaz * irho
            
            Lw(i,j,5,1) = 0.d0
            Lw(i,j,5,2) = 0.d0
            Lw(i,j,5,3) = 0.5d0 * betaz
            Lw(i,j,5,4) = -0.5d0 * betay
            Lw(i,j,5,5) = -0.5d0 * betaz * S * irhos
            Lw(i,j,5,6) = 0.5d0 * betay * S * irhos
            
            Lw(i,j,6,1) = Nf * alphaf * asqr * irho
            Lw(i,j,6,2) = Nf * Cff
            Lw(i,j,6,3) = -Nf * Qs * betay
            Lw(i,j,6,4) = -Nf * Qs * betaz
            Lw(i,j,6,5) = Nf * As * betay * irho
            Lw(i,j,6,6) = Nf * As * betaz * irho
            
        END DO
    END DO
    
    

END SUBROUTINE  compute_primitive_eigenvectors


SUBROUTINE compute_interface_states_prim(ilow, ihi, jlow, jhi)

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi
    INTEGER :: i, j, k, offset, s_i(1-nb:nmax+nb,1-nb:nmax+nb)
    REAL*8 :: delw_L(6), delw_R(6), delw_C(6), delw_0(6), delw_lim(6), w_R(6), wL_tmp, wR_tmp, mabs_slope, &
              delw(1-nb:nmax+nb,1-nb:nmax+nb,6), wL_lin(1-nb:nmax+nb,1-nb:nmax+nb,6), wR_lin(1-nb:nmax+nb,1-nb:nmax+nb,6)
    REAL*8 :: signL, signR, signC, third, sixth, s112, s712,  tmp1, tmp2, tmp3, eta(6), &
              eta_tilde(6), d2w(6), d2w_L(6), d2w_R(6), flat_coeff, w_sh, drho, dv, drho2, drho_rat, &
              flat_tilde(1-nb:nmax+nb,1-nb:nmax+nb)
    REAL*8, PARAMETER :: eta1 = 20.d0 , eta2 = 0.005d0 , eps = 0.01, omega1 = 0.75d0, omega2 = 10.d0
    
    third = 1.d0/3.d0
    sixth = 1.d0/6.d0
    s112 = 1.d0/12.d0
    s712 = 7.d0/12.d0
   
    IF(.NOT. piecewise_parabolic) THEN
    
        DO j = jlow-nb, jhi+nb
            DO i = ilow-nb+1, ihi+nb-1

                ! Compute first order (cell-center) left, right and centered differences of the primitive variables 
                delw_L(:) = w(i,j,:) - w(i-1,j,:) 
                delw_R(:) = w(i+1,j,:) - w(i,j,:) 
                delw_C(:) = 0.5d0 *( delw_L(:) + delw_R(:) )

                ! compute monotonized slope of primitive variables
                DO k = 1, 6
                    
                    ! MC (monotonized central-difference) limiter
                    delw(i,j,k) = SIGN(1.d0, delw_C(k)) * MIN( MC_beta * ABS(delw_L(k)), MC_beta * ABS(delw_R(k)), ABS(delw_C(k)) )                
                        
                    ! min-mod limiter (most stable)
                    !delw(i,j,k) =0.5d0 * (SIGN(1.d0, delw_L(k)) + SIGN(1.d0, delw_R(k))) * MIN(ABS(delw_L(k)), ABS(delw_R(k)))                        
                
                END DO
                
            END DO
        END DO

        ! compute left and right interface states
        DO j = jlow-nb, jhi+nb
            DO i = ilow-nb+1, ihi+nb-2

                ! left state
                qintL(i,j,1) =  w(i,j,1) + 0.5d0 * delw(i,j,1)  
                qintL(i,j,2) = (w(i,j,2) + 0.5d0 * delw(i,j,2)) ! * qintL(i,j,1)  
                qintL(i,j,3) = (w(i,j,3) + 0.5d0 * delw(i,j,3)) ! * qintL(i,j,1) 
                qintL(i,j,4) = (w(i,j,4) + 0.5d0 * delw(i,j,4)) ! * qintL(i,j,1) 
                qintL(i,j,5) =  wb(i,j)   
                qintL(i,j,6) =  w(i,j,5) + 0.5d0 * delw(i,j,5)  
                qintL(i,j,7) =  w(i,j,6) + 0.5d0 * delw(i,j,6)  

            END DO
        END DO

        DO j = jlow-nb, jhi+nb
            DO i = ilow-nb+1, ihi+nb-2

                ! right state
                qintR(i,j,1) =  w(i+1,j,1) - 0.5d0 * delw(i+1,j,1)  
                qintR(i,j,2) = (w(i+1,j,2) - 0.5d0 * delw(i+1,j,2)) ! * qintR(i,j,1)  
                qintR(i,j,3) = (w(i+1,j,3) - 0.5d0 * delw(i+1,j,3)) ! * qintR(i,j,1) 
                qintR(i,j,4) = (w(i+1,j,4) - 0.5d0 * delw(i+1,j,4)) ! * qintR(i,j,1) 
                qintR(i,j,5) =  wb(i,j)  
                qintR(i,j,6) =  w(i+1,j,5) - 0.5d0 * delw(i+1,j,5)  
                qintR(i,j,7) =  w(i+1,j,6) - 0.5d0 * delw(i+1,j,6)  

            END DO
        END DO
    
        offset = 0
        
    ELSE

        DO j = jlow-nb, jhi+nb
            DO i = ilow-nb+1, ihi+nb-1

                ! Compute 2nd order (cell-center) left, right and centered differences of the primitive variables 
                delw_L(:) = w(i,j,:) - w(i-1,j,:) 
                delw_R(:) = w(i+1,j,:) -w(i,j,:)  
                delw_C(:) = 0.5d0 *(delw_L(:) + delw_R(:))
                
                delw_lim(:) = 0.d0
                
                ! compute monotonized slope of primitive variables using MC (monotonized central-difference) beta limiter
                ! (this step will ensure that the cell edge value lies between the neighboring cell celter values, i.e. no new extrema will be created)
                DO k = 1, 6
                    
                    mabs_slope = MIN(ABS(delw_L(k)), ABS(delw_R(k)))
                    
                    ! check whether left and right differences have the same sign 
                    IF((delw_L(k)*delw_R(k)) .GT. 0.d0) THEN
                        
                        delw_lim(k) = SIGN(1.d0, delw_C(k)) * MIN(MC_beta*mabs_slope, ABS(delw_C(k)))                
                        !delw_lim(k) = 0.5d0 * (SIGN(1.d0, delw_L(k)) + SIGN(1.d0, delw_R(k))) * MIN(MC_beta*mabs_slope, ABS(delw_C(k)))                
                    
                    END IF                     
                    
                END DO
                
                delw(i,j,:) = delw_lim(:)                                 
                    
            END DO
        END DO

        ! quadratic interpolation of left and right interface values 
        DO j = jlow-nb, jhi+nb
            DO i = ilow-nb+2, ihi+nb-2

                wL(i,j,:) = 0.5d0 * (w(i-1,j,:) + w(i,j,:)) + sixth * (delw(i-1,j,:)-delw(i,j,:)) 
                wR(i,j,:) = 0.5d0 * (w(i,j,:) + w(i+1,j,:)) + sixth * (delw(i,j,:)-delw(i+1,j,:)) 
                 
            END DO
        END DO
        
        !******************************************************************************************************
        ! Discontinuiuty Steepener (see WC eqn. 1.15) 
        !******************************************************************************************************
        IF(discontinuiuty_steepener) THEN        
      
          ! linear interpolation of left and right interface values
            DO j = jlow-nb, jhi+nb
                DO i = ilow-nb+2, ihi+nb-2

                    wL_lin(i,j,:) = w(i-1,j,:) + 0.5d0 * delw(i-1,j,:)
                    wR_lin(i,j,:) = w(i+1,j,:) - 0.5d0 * delw(i+1,j,:) 
                     
                END DO
            END DO

            
            ! Detect discontinuiuties and semoothly switch from piecewise quadratic interpolated interface states
            ! to a piecewise linear interpolanted values. (Piecewise linear should yield a sharper profile near discontinuities) 
            DO j = jlow-nb, jhi+nb
                DO i = ilow-nb+2, ihi+nb-2

                    ! approximate the second derivative at the cell edges
                    d2w_L(:) = sixth * ((w(i+2,j,:)-w(i+1,j,:)) - (w(i+1,j,:)-w(i,j,:)))
                    d2w_R(:) = sixth * ((w(i,j,:)-w(i-1,j,:)) - (w(i-1,j,:)-w(i-2,j,:)))
                    d2w(:) = w(i+1,j,:) - w(i-1,j,:)
                    
                    eta_tilde(:) = 0.d0
     
                    ! detect discontinuiuty and compute piecewise linear weighting factor ("eta")
                    DO k = 1, 6 
                    
                        tmp1 = eps * MIN(ABS(w(i+1,j,k)), ABS(w(i-1,j,k)))
                    
                        ! detect a discontinuiuty by checking if second derivative changes sign and third derivative becomes "large" (largeness determined by the eps parameter)
                        IF( ((d2w_L(k)*d2w_R(k)) .LT. 0.d0) .AND. ((ABS(d2w(k))-tmp1) .GT. 0.d0)) THEN
                            eta_tilde(k) =  -(d2w_R(k)-d2w_L(k)) / d2w(k)
                        END IF
                   
                        eta(k) = steepen_var(k) * MAX( 0.d0, MIN( eta1*(eta_tilde(k)-eta2), 1.d0) )
                   
                    END DO
                    
                    ! compute modified left and right interface states (weighted average of linear and quadratic interpolants)
                    wL(i,j,:) = (1.d0 - eta(:)) * wL(i,j,:) + eta(:) * wL_lin(i,j,:) 
                    wR(i,j,:) = (1.d0 - eta(:)) * wR(i,j,:) + eta(:) * wR_lin(i,j,:) 
                    
                END DO
            END DO

        END IF

        
        ! apply further monotonicity constraints to ensure that w_i between wL_i and wR_i 
        DO j = jlow-nb, jhi+nb
            DO i = ilow-nb+2, ihi+nb-2
                DO k = 1, 6
                
                    tmp1 = (wR(i,j,k)-w(i,j,k)) * (w(i,j,k)-wL(i,j,k))
                    tmp2 = 6.d0 * (wR(i,j,k)-wL(i,j,k)) * (w(i,j,k)-0.5d0*(wL(i,j,k)+wR(i,j,k)))
                    tmp3 = (wR(i,j,k)-wL(i,j,k))**2

                    IF(tmp1 .LE. 0.d0) THEN  ! if we're at a local extremum, then we completely flatten the parabola (so that left and right edge values are equal)  
                        wL(i,j,k) = w(i,j,k) - 0.5d0 * delw(i,j,k)
                        wR(i,j,k) = w(i,j,k) + 0.5d0 * delw(i,j,k)
                    ELSE IF(tmp2 .GT. tmp3) THEN ! reduce curvature of the parabola
                        wL(i,j,k) = 3.d0*w(i,j,k) - 2.d0*wR(i,j,k)
                    ELSE IF(tmp2 .LT. -tmp3) THEN
                        wR(i,j,k) = 3.d0*w(i,j,k) - 2.d0*wL(i,j,k)
                    END IF
                    
                END DO
            END DO
        END DO        
        
        
        !************************************************************************************************************  
        
        GO TO 1000
        
        ! compute improved estimates for monotonized slopes and interface states
        DO j = jlow-nb, jhi+nb
            DO i = ilow-nb+2, ihi+nb-2
            
                delw_lim(:) = 0.d0

                DO k = 1, 6

                    tmp1 = (wR(i,j,k)-w(i,j,k)) * (w(i,j,k)-wL(i,j,k))
                    tmp2 = (w(i-1,j,k)-w(i,j,k)) * (w(i,j,k)-w(i+1,j,k))
                    
                    IF((tmp1 .LE. 0.d0) .AND. (tmp2 .LE. 0.d0)) THEN
                    
                        ! Compute 2nd order (cell-center) left, right and centered differences of the primitive variables 
                        delw_L(k) = w(i-2,j,k) - 2.d0 * w(i-1,j,k) + w(i,j,k)   
                        delw_R(k) = w(i,j,k) - 2.d0 * w(i+1,j,k)   + w(i+2,j,k)  
                        delw_C(k) = w(i-1,j,k) - 2.d0 * w(i,j,k)   + w(i+1,j,k)
                        delw_0(k) = -12.d0 * (w(i,j,k) - 0.5d0 * (wL(i,j,k) + wR(i,j,k)))
                                                
                         ! compute monotonized slope of primitive variables using MC (monotonized central-difference) beta limiter    
                        mabs_slope = MIN(ABS(delw_L(k)), ABS(delw_R(k)), ABS(delw_C(k)))
                            
                        ! check whether all four differences have the same sign 
                        IF( ((delw_L(k) .GT. 0.d0) .AND. (delw_C(k) .GT. 0.d0) .AND. (delw_R(k) .GT. 0.d0) .AND. (delw_0(k) .GT. 0.d0)) .OR. &
                            ((delw_L(k) .LT. 0.d0) .AND. (delw_C(k) .LT. 0.d0) .AND. (delw_R(k) .LT. 0.d0) .AND. (delw_0(k) .LT. 0.d0))) THEN
                                
                            delw_lim(k) = SIGN(1.d0, delw_0(k)) * MIN(MC_beta*mabs_slope, ABS(delw_0(k)))                
                           
                        END IF

                        IF(ABS(delw_0(k)) .GT. 0.d0) THEN
                            delw_lim(k) = delw_lim(k) / delw_0(k) 
                        ELSE
                            delw_lim(k) = 0.d0
                        END IF
                   
                        wL(i,j,k) =  w(i,j,k) + (wL(i,j,k) - w(i,j,k)) * delw_lim(k)
                        wR(i,j,k) =  w(i,j,k) + (wR(i,j,k) - w(i,j,k)) * delw_lim(k)        
                                      
                   END IF                     
 
                END DO        
                                
            END DO
        END DO


       ! apply further monotonicity constraint to ensure these values lie between neighboring cell center values (NOTE: TRY PUTTING THIS AFTER THE NEXT LOOP AND REPLACE wL/R with qintL/R)
        DO j = jlow-nb, jhi+nb
            DO i = ilow-nb+2, ihi+nb-2
                DO k = 1, 6

                    tmp1 = (wR(i,j,k)-w(i,j,k)) * (w(i,j,k)-wL(i,j,k))
                    tmp2 = 6.d0 * (wR(i,j,k)-wL(i,j,k)) * (w(i,j,k)-0.5d0*(wL(i,j,k)+wR(i,j,k)))
                    tmp3 = (wR(i,j,k)-wL(i,j,k))**2

                    IF(tmp1 .LE. 0.d0) THEN   
                        wL(i,j,k) = w(i,j,k)
                        wR(i,j,k) = w(i,j,k)
                    ELSE IF(tmp2 .GT. tmp3) THEN
                        wL(i,j,k) = 3.d0*w(i,j,k) - 2.d0*wR(i,j,k)
                    ELSE IF(tmp2 .LT. -tmp3) THEN
                        wR(i,j,k) = 3.d0*w(i,j,k) - 2.d0*wL(i,j,k)
                    END IF
                                
                END DO    
            END DO
        END DO

        1000 CONTINUE
     

        ! compute left and right interface states
        DO j = jlow-nb, jhi+nb
            DO i = ilow-nb+2, ihi+nb-3

                ! left state
                qintL(i,j,1) = wR(i,j,1) 
                qintL(i,j,2) = wR(i,j,2) ! * qintL(i,j,1) 
                qintL(i,j,3) = wR(i,j,3) ! * qintL(i,j,1) 
                qintL(i,j,4) = wR(i,j,4) ! * qintL(i,j,1) 
                qintL(i,j,5) = wb(i,j)   
                qintL(i,j,6) = wR(i,j,5)  
                qintL(i,j,7) = wR(i,j,6)  

            END DO
        END DO

        DO j = jlow-nb, jhi+nb
            DO i = ilow-nb+2, ihi+nb-3

                ! right state
                qintR(i,j,1) = wL(i+1,j,1)  
                qintR(i,j,2) = wL(i+1,j,2) ! * qintR(i,j,1)  
                qintR(i,j,3) = wL(i+1,j,3) ! * qintR(i,j,1) 
                qintR(i,j,4) = wL(i+1,j,4) ! * qintR(i,j,1) 
                qintR(i,j,5) = wb(i,j)  
                qintR(i,j,6) = wL(i+1,j,5)  
                qintR(i,j,7) = wL(i+1,j,6)  

            END DO
        END DO
        
        offset = 1

    END IF

   
    !GO TO 1001
    ! explicitly check for and enforce monotonicity, i.e. make sure that MIN(w_i,w_i+1) <= (wL_i+1/2, wR_i+1/2) <= MAX(w_i,w_i+1)
    !DO j = jlow-nb, jhi+nb
    !    DO i = ilow-nb+1+offset, ihi+nb-2-offset
    !
    !        qintL(i,j,1:4) = MAX(MIN(w(i,j,1:4), w(i+1,j,1:4)), qintL(i,j,1:4))
    !        qintL(i,j,1:4) = MIN(MAX(w(i,j,1:4), w(i+1,j,1:4)), qintL(i,j,1:4))
    !        qintL(i,j,6:7) = MAX(MIN(w(i,j,5:6), w(i+1,j,5:6)), qintL(i,j,6:7))
    !        qintL(i,j,6:7) = MIN(MAX(w(i,j,5:6), w(i+1,j,5:6)), qintL(i,j,6:7))
    !        
    !        qintR(i,j,1:4) = MAX(MIN(w(i,j,1:4), w(i+1,j,1:4)), qintR(i,j,1:4))
    !        qintR(i,j,1:4) = MIN(MAX(w(i,j,1:4), w(i+1,j,1:4)), qintR(i,j,1:4))            
    !        qintR(i,j,6:7) = MAX(MIN(w(i,j,5:6), w(i+1,j,5:6)), qintR(i,j,6:7))
    !        qintR(i,j,6:7) = MIN(MAX(w(i,j,5:6), w(i+1,j,5:6)), qintR(i,j,6:7))
    !        
    !    END DO
    !END DO
    !1001 CONTINUE    
        
    ! convert from primitive to conserved variables
    DO j = jlow-nb, jhi+nb
        DO i = ilow-nb+1+offset, ihi+nb-2-offset

            ! left state
            qintL(i,j,2) = qintL(i,j,2) * qintL(i,j,1)  
            qintL(i,j,3) = qintL(i,j,3) * qintL(i,j,1) 
            qintL(i,j,4) = qintL(i,j,4) * qintL(i,j,1) 

            ! right state
            qintR(i,j,2) = qintR(i,j,2) * qintR(i,j,1)  
            qintR(i,j,3) = qintR(i,j,3) * qintR(i,j,1) 
            qintR(i,j,4) = qintR(i,j,4) * qintR(i,j,1) 

        END DO
    END DO    
        
    
    
END SUBROUTINE compute_interface_states_prim


SUBROUTINE compute_interface_states_char(ilow, ihi, jlow, jhi)

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi
    INTEGER :: i, j, k, offset
    REAL*8 :: delw_L(6), delw_R(6), delw_C(6), delw_0(6), delw_lim(6), w_R(6), wL_tmp, wR_tmp, mabs_slope, &
              delw(1-nb:nmax+nb,1-nb:nmax+nb,6), wL_lin(1-nb:nmax+nb,1-nb:nmax+nb,6), wR_lin(1-nb:nmax+nb,1-nb:nmax+nb,6)
    REAL*8 :: signL, signR, signC, third, sixth, s112, s712,  tmp1, tmp2, tmp3, eta(6), &
              eta_tilde(6), d2w(6), d2w_L(6), d2w_R(6)
    REAL*8, PARAMETER :: eta1 = 20.d0 , eta2 = 0.005d0 , eps = 0.01, omega1 = 0.75d0, omega2 = 10.d0
    
    REAL*8 :: a_C(6), a_L_1(6), a_R_1(6), a_R_2(6), dela_L_1(6), dela_L_2(6), dela_R_1(6), dela_C_1(6), &
              dela_R_2(6), dela_C_2(6), dela_m(6), dela_m_1(6), dela_m_2(6)
    
    
    third = 1.d0/3.d0
    sixth = 1.d0/6.d0
    s112 = 1.d0/12.d0
    s712 = 7.d0/12.d0
   
    IF(.NOT. piecewise_parabolic) THEN
    
        DO j = jlow-nb, jhi+nb
            DO i = ilow-nb+1, ihi+nb-1


               ! locally project the characteristic variables
                DO k = 1, 6     

                    a_L_1(k) = Lw(i,j,k,1) * w(i-1,j,1) + Lw(i,j,k,2) * w(i-1,j,2) + Lw(i,j,k,3) * w(i-1,j,3) + & 
                               Lw(i,j,k,4) * w(i-1,j,4) + Lw(i,j,k,5) * w(i-1,j,5) + Lw(i,j,k,6) * w(i-1,j,6)
                               
                    a_C(k)   = Lw(i,j,k,1) * w(i,j,1) + Lw(i,j,k,2) * w(i,j,2) + Lw(i,j,k,3) * w(i,j,3) + & 
                               Lw(i,j,k,4) * w(i,j,4) + Lw(i,j,k,5) * w(i,j,5) + Lw(i,j,k,6) * w(i,j,6)
                   
                    a_R_1(k) = Lw(i,j,k,1) * w(i+1,j,1) + Lw(i,j,k,2) * w(i+1,j,2) + Lw(i,j,k,3) * w(i+1,j,3) + & 
                               Lw(i,j,k,4) * w(i+1,j,4) + Lw(i,j,k,5) * w(i+1,j,5) + Lw(i,j,k,6) * w(i+1,j,6)
             
                END DO
            
                ! Compute (cell-center) left, right and centered differences of the characteristic variables 
                dela_L_1(:) = a_C(:) - a_L_1(:) 
                dela_R_1(:) = a_R_1(:) - a_C(:)  
                dela_C_1(:) = 0.5d0 *( dela_L_1(:) + dela_R_1(:) )

                dela_m_1(:) = 0.d0

                ! compute monotonized slope of characteristic variables
                DO k = 1, 6
                    
                    ! check whether left and right differences have the same sign and apply limiter
                    IF((dela_L_1(k)*dela_R_1(k)) .GT. 0.d0) THEN
                        
                        dela_m_1(k) = SIGN(1.d0, dela_C_1(k)) * MIN(MC_beta*ABS(dela_L_1(k)), MC_beta*ABS(dela_R_1(k)), ABS(dela_C_1(k)))                
                    
                    END IF        
                    
                       
                END DO               
                
                ! project the monotonized slope back onto primitive variables
                DO k = 1, 6
                
                    delw(i,j,k)  = dela_m_1(1) * Rw(i,j,k,1) + dela_m_1(2) * Rw(i,j,k,2) + dela_m_1(3) * Rw(i,j,k,3) + & 
                                   dela_m_1(4) * Rw(i,j,k,4) + dela_m_1(5) * Rw(i,j,k,5) + dela_m_1(6) * Rw(i,j,k,6)                        
                                   
                END DO
              

            END DO
        END DO

        ! compute left and right interface states
        DO j = jlow-nb, jhi+nb
            DO i = ilow-nb+1, ihi+nb-2

                ! left state
                qintL(i,j,1) =  w(i,j,1) + 0.5d0 * delw(i,j,1)  
                qintL(i,j,2) = (w(i,j,2) + 0.5d0 * delw(i,j,2)) 
                qintL(i,j,3) = (w(i,j,3) + 0.5d0 * delw(i,j,3)) 
                qintL(i,j,4) = (w(i,j,4) + 0.5d0 * delw(i,j,4)) 
                qintL(i,j,5) =  wb(i,j)   
                qintL(i,j,6) =  w(i,j,5) + 0.5d0 * delw(i,j,5)  
                qintL(i,j,7) =  w(i,j,6) + 0.5d0 * delw(i,j,6)  

            END DO
        END DO

        DO j = jlow-nb, jhi+nb
            DO i = ilow-nb+1, ihi+nb-2

                ! right state
                qintR(i,j,1) =  w(i+1,j,1) - 0.5d0 * delw(i+1,j,1)  
                qintR(i,j,2) = (w(i+1,j,2) - 0.5d0 * delw(i+1,j,2))   
                qintR(i,j,3) = (w(i+1,j,3) - 0.5d0 * delw(i+1,j,3)) 
                qintR(i,j,4) = (w(i+1,j,4) - 0.5d0 * delw(i+1,j,4))  
                qintR(i,j,5) =  wb(i,j)  
                qintR(i,j,6) =  w(i+1,j,5) - 0.5d0 * delw(i+1,j,5)  
                qintR(i,j,7) =  w(i+1,j,6) - 0.5d0 * delw(i+1,j,6)  

            END DO
        END DO
    
        offset = 0
        
   
   ELSE

        DO j = jlow-nb, jhi+nb
            DO i = ilow-nb+1, ihi+nb-2
            
                ! locally project the characteristic variables
                DO k = 1, 6     

                    a_L_1(k) = Lw(i,j,k,1) * w(i-1,j,1) + Lw(i,j,k,2) * w(i-1,j,2) + Lw(i,j,k,3) * w(i-1,j,3) + & 
                               Lw(i,j,k,4) * w(i-1,j,4) + Lw(i,j,k,5) * w(i-1,j,5) + Lw(i,j,k,6) * w(i-1,j,6)
                               
                    a_C(k)   = Lw(i,j,k,1) * w(i,j,1) + Lw(i,j,k,2) * w(i,j,2) + Lw(i,j,k,3) * w(i,j,3) + & 
                               Lw(i,j,k,4) * w(i,j,4) + Lw(i,j,k,5) * w(i,j,5) + Lw(i,j,k,6) * w(i,j,6)
                   
                    a_R_1(k) = Lw(i,j,k,1) * w(i+1,j,1) + Lw(i,j,k,2) * w(i+1,j,2) + Lw(i,j,k,3) * w(i+1,j,3) + & 
                               Lw(i,j,k,4) * w(i+1,j,4) + Lw(i,j,k,5) * w(i+1,j,5) + Lw(i,j,k,6) * w(i+1,j,6)

                    a_R_2(k) = Lw(i,j,k,1) * w(i+2,j,1) + Lw(i,j,k,2) * w(i+2,j,2) + Lw(i,j,k,3) * w(i+2,j,3) + & 
                               Lw(i,j,k,4) * w(i+2,j,4) + Lw(i,j,k,5) * w(i+2,j,5) + Lw(i,j,k,6) * w(i+2,j,6)
                               
                END DO

            
                ! Compute (cell-center) left, right and centered differences of the characteristic variables 
                dela_L_1(:) = a_C(:) - a_L_1(:) 
                dela_R_1(:) = a_R_1(:) - a_C(:)  
                dela_C_1(:) = 0.5d0 *( dela_L_1(:) + dela_R_1(:) )

                dela_L_2(:) = a_R_1(:) - a_C(:) 
                dela_R_2(:) = a_R_2(:) - a_R_1(:)  
                dela_C_2(:) = 0.5d0 *( dela_L_2(:) + dela_R_2(:) )

                dela_m_1(:) = 0.d0
                dela_m_2(:) = 0.d0

                ! compute monotonized slope of characteristic variables
                DO k = 1, 6
                    
                   ! check whether left and right differences have the same sign and apply limiter
                    IF((dela_L_1(k)*dela_R_1(k)) .GT. 0.d0) THEN
                        
                        dela_m_1(k) = SIGN(1.d0, dela_C_1(k)) * MIN(MC_beta*ABS(dela_L_1(k)), MC_beta*ABS(dela_R_1(k)), ABS(dela_C_1(k)))                
                    
                    END IF        
                    
                    ! check whether left and right differences have the same sign and apply limiter
                    IF((dela_L_2(k)*dela_R_2(k)) .GT. 0.d0) THEN
                        
                        dela_m_2(k) = SIGN(1.d0, dela_C_2(k)) * MIN(MC_beta*ABS(dela_L_2(k)), MC_beta*ABS(dela_R_2(k)), ABS(dela_C_2(k)))                
                    
                    END IF  

                       
                END DO
                
                dela_m(:) = dela_m_2(:) - dela_m_1(:) 
                
                ! project the monotonized slope back onto primitive variables
                DO k = 1, 6
                
                    delw(i,j,k)  = dela_m(1) * Rw(i,j,k,1) + dela_m(2) * Rw(i,j,k,2) + dela_m(3) * Rw(i,j,k,3) + & 
                                   dela_m(4) * Rw(i,j,k,4) + dela_m(5) * Rw(i,j,k,5) + dela_m(6) * Rw(i,j,k,6)                    
                                   
                END DO
                    
            END DO
        END DO

        ! quadratic interpolation of left and right interface values 
        DO j = jlow-nb, jhi+nb
            DO i = ilow-nb+2, ihi+nb-2

                wL(i,j,:) = 0.5d0 * (w(i-1,j,:) + w(i,j,:)) - sixth * delw(i-1,j,:) 
                wR(i,j,:) = 0.5d0 * (w(i,j,:) + w(i+1,j,:)) - sixth * delw(i,j,:) 
                 
            END DO
        END DO
   
        !******************************************************************************************************
        ! Discontinuiuty Steepener (see WC eqn. 1.15) 
        !******************************************************************************************************
        IF(discontinuiuty_steepener) THEN        

            ! linear interpolation of left and right interface values
            DO j = jlow-nb, jhi+nb
                DO i = ilow-nb+2, ihi+nb-2

                    wL_lin(i,j,:) = w(i-1,j,:) + 0.5d0 * delw(i-1,j,:)
                    wR_lin(i,j,:) = w(i+1,j,:) - 0.5d0 * delw(i+1,j,:) 
                     
                END DO
            END DO

            
            ! Detect discontinuiuties and semoothly switch from piecewise quadratic interpolated interface states
            ! to a piecewise linear interpolanted values. (Piecewise linear should yield a sharper profile near discontinuities) 
            DO j = jlow-nb, jhi+nb
                DO i = ilow-nb+2, ihi+nb-2

                    ! approximate the second derivative at the cell edges
                    d2w_L(:) = sixth * ((w(i+2,j,:)-w(i+1,j,:)) - (w(i+1,j,:)-w(i,j,:)))
                    d2w_R(:) = sixth * ((w(i,j,:)-w(i-1,j,:)) - (w(i-1,j,:)-w(i-2,j,:)))
                    d2w(:) = w(i+1,j,:) - w(i-1,j,:)
                    
                    eta_tilde(:) = 0.d0
                    eta(:) = 0.d0
     
                    ! detect discontinuiuty and compute piecewise linear weighting factor ("eta")
                    DO k = 1, 6 
                    
                        tmp1 = eps * MIN(ABS(w(i+1,j,k)), ABS(w(i-1,j,k)))
                    
                        ! detect a discontinuiuty by checking if second derivative changes sign and third derivative becomes "large" (largeness determined by the eps parameter)
                        IF( ((d2w_L(k)*d2w_R(k)) .LT. 0.d0) .AND. ((ABS(d2w(k))-tmp1) .GT. 0.d0)) THEN
                            eta_tilde(k) =  -(d2w_R(k)-d2w_L(k)) / d2w(k)
                        END IF
                   
                        eta(k) = steepen_var(k) * MAX( 0.d0, MIN( eta1*(eta_tilde(k)-eta2), 1.d0) )
                   
                    END DO
                    
                    ! compute modified left and right interface states (weighted average of linear and quadratic interpolants)
                    wL(i,j,:) = (1.d0 - eta(:)) * wL(i,j,:) + eta(:) * wL_lin(i,j,:) 
                    wR(i,j,:) = (1.d0 - eta(:)) * wR(i,j,:) + eta(:) * wR_lin(i,j,:) 
                    
                END DO
            END DO

        END IF
        
        ! apply further monotonicity constraints to ensure that w_i between wL_i and wR_i 
        DO j = jlow-nb, jhi+nb
            DO i = ilow-nb+2, ihi+nb-2
                DO k = 1, 6
                
                    tmp1 = (wR(i,j,k)-w(i,j,k)) * (w(i,j,k)-wL(i,j,k))
                    tmp2 = 6.d0 * (wR(i,j,k)-wL(i,j,k)) * (w(i,j,k)-0.5d0*(wL(i,j,k)+wR(i,j,k)))
                    tmp3 = (wR(i,j,k)-wL(i,j,k))**2

                    IF(tmp1 .LE. 0.d0) THEN  ! if we're at a local extremum, then we completely flatten the parabola (so that left and right edge values are equal)  
                        wL(i,j,k) = w(i,j,k) 
                        wR(i,j,k) = w(i,j,k) 
                    ELSE IF(tmp2 .GT. tmp3) THEN ! this part here will reduce the curvature of the parabola
                        wL(i,j,k) = 3.d0*w(i,j,k) - 2.d0*wR(i,j,k)
                    ELSE IF(tmp2 .LT. -tmp3) THEN
                        wR(i,j,k) = 3.d0*w(i,j,k) - 2.d0*wL(i,j,k)
                    END IF
                    
                END DO
            END DO
        END DO        
        
        ! compute left and right interface states
        DO j = jlow-nb, jhi+nb
            DO i = ilow-nb+2, ihi+nb-3

                ! left state
                qintL(i,j,1) = wR(i,j,1) 
                qintL(i,j,2) = wR(i,j,2) 
                qintL(i,j,3) = wR(i,j,3)  
                qintL(i,j,4) = wR(i,j,4)  
                qintL(i,j,5) = wb(i,j)   
                qintL(i,j,6) = wR(i,j,5)  
                qintL(i,j,7) = wR(i,j,6)  

            END DO
        END DO

        DO j = jlow-nb, jhi+nb
            DO i = ilow-nb+2, ihi+nb-3

                ! right state
                qintR(i,j,1) = wL(i+1,j,1)  
                qintR(i,j,2) = wL(i+1,j,2)  
                qintR(i,j,3) = wL(i+1,j,3) 
                qintR(i,j,4) = wL(i+1,j,4)  
                qintR(i,j,5) = wb(i,j)  
                qintR(i,j,6) = wL(i+1,j,5)  
                qintR(i,j,7) = wL(i+1,j,6)  

            END DO
        END DO
        
        offset = 1

    END IF

   
    ! explicitly check for and enforce monotonicity, i.e. make sure that MIN(w_i,w_i+1) <= (wL_i+1/2, wR_i+1/2) <= MAX(w_i,w_i+1)
    !DO j = jlow-nb, jhi+nb
    !    DO i = ilow-nb+1+offset, ihi+nb-2-offset
    !
    !        qintL(i,j,1:4) = MAX(MIN(w(i,j,1:4), w(i+1,j,1:4)), qintL(i,j,1:4))
    !        qintL(i,j,1:4) = MIN(MAX(w(i,j,1:4), w(i+1,j,1:4)), qintL(i,j,1:4))
    !        qintL(i,j,6:7) = MAX(MIN(w(i,j,5:6), w(i+1,j,5:6)), qintL(i,j,6:7))
    !        qintL(i,j,6:7) = MIN(MAX(w(i,j,5:6), w(i+1,j,5:6)), qintL(i,j,6:7))
    !        
    !        qintR(i,j,1:4) = MAX(MIN(w(i,j,1:4), w(i+1,j,1:4)), qintR(i,j,1:4))
    !        qintR(i,j,1:4) = MIN(MAX(w(i,j,1:4), w(i+1,j,1:4)), qintR(i,j,1:4))            
    !        qintR(i,j,6:7) = MAX(MIN(w(i,j,5:6), w(i+1,j,5:6)), qintR(i,j,6:7))
    !        qintR(i,j,6:7) = MIN(MAX(w(i,j,5:6), w(i+1,j,5:6)), qintR(i,j,6:7))
    !        
    !    END DO
    !END DO

        
    ! convert from primitive to conserved variables
    DO j = jlow-nb, jhi+nb
        DO i = ilow-nb+1+offset, ihi+nb-2-offset

            ! left state
            qintL(i,j,2) = qintL(i,j,2) * qintL(i,j,1)  
            qintL(i,j,3) = qintL(i,j,3) * qintL(i,j,1) 
            qintL(i,j,4) = qintL(i,j,4) * qintL(i,j,1) 

            ! right state
            qintR(i,j,2) = qintR(i,j,2) * qintR(i,j,1)  
            qintR(i,j,3) = qintR(i,j,3) * qintR(i,j,1) 
            qintR(i,j,4) = qintR(i,j,4) * qintR(i,j,1) 

        END DO
    END DO    
        
    
    
END SUBROUTINE compute_interface_states_char



! Compute the cell-centered refernece emf
SUBROUTINE compute_cell_center_emf_3d(ilow, ihi, jlow, jhi, klow, khi) 

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi, klow, khi
    INTEGER :: i, j, k
    REAL*8 :: irho
    
           
    DO k = klow-nb, khi+nb
        DO j = jlow-nb, jhi+nb
            DO i = ilow-nb, ihi+nb         
       
                irho = 1.d0 / q_3d(i,j,k,1)
                emf_3d(i,j,k,1,3) = ( q_3d(i,j,k,6) * q_3d(i,j,k,4)  - q_3d(i,j,k,7) * q_3d(i,j,k,3) ) * irho   ! Ex_i,j,k          
                emf_3d(i,j,k,2,3) = ( q_3d(i,j,k,7) * q_3d(i,j,k,2)  - q_3d(i,j,k,5) * q_3d(i,j,k,4) ) * irho   ! Ey_i,j,k      
                emf_3d(i,j,k,3,3) = ( q_3d(i,j,k,5) * q_3d(i,j,k,3)  - q_3d(i,j,k,6) * q_3d(i,j,k,2) ) * irho   ! Ez_i,j,k 
                
            END DO
        END DO    
    END DO


    IF(flux_protection) THEN
        
        DO k = klow-nb, khi+nb
            DO j = jlow-nb, jhi+nb
                DO i = ilow-nb, ihi+nb     

                    emf_3d_prot(i,j,k,1,3) = emf_3d(i,j,k,1,3)
                    emf_3d_prot(i,j,k,2,3) = emf_3d(i,j,k,2,3)
                    emf_3d_prot(i,j,k,3,3) = emf_3d(i,j,k,3,3)
                        

                END DO
            END DO    
        END DO
        
    END IF
    

END SUBROUTINE compute_cell_center_emf_3d



! Compute the cell-corner emf

!**********************************************************************
! Electric field array layout: emf_3d(i,j,k,1,1)   = Ex_i,j+1/2,k
!                              emf_3d(i,j,k,1,2)   = Ex_i,j,k+1/2  
!                              emf_3d(i,j,k,1,3)   = Ex_i,j,k  
!                              emf_corner(i,j,k,1) = Ex_i,j+1/2,k+1/2
!
!                              emf_3d(i,j,k,2,1)   = Ey_i+1/2,j,k
!                              emf_3d(i,j,k,2,2)   = Ey_i,j,k+1/2  
!                              emf_3d(i,j,k,2,3)   = Ey_i,j,k  
!                              emf_corner(i,j,k,2) = Ey_i+1/2,j,k+1/2
!
!                              emf_3d(i,j,k,3,1)   = Ez_i+1/2,j,k
!                              emf_3d(i,j,k,3,2)   = Ez_i,j+1/2,k  
!                              emf_3d(i,j,k,3,3)   = Ez_i,j,k 
!                              emf_corner(i,j,k,3) = Ez_i+1/2,j+1/2,k 
!**********************************************************************
SUBROUTINE compute_cell_corner_emf_3d(emf_3d_in, emf_corner_in, ilow, ihi, jlow, jhi, klow, khi)

    REAL*8, INTENT(IN) :: emf_3d_in(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,3,3)
    REAL*8, INTENT(INOUT) :: emf_corner_in(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,3)
    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi, klow, khi
    INTEGER :: i, j, k
    REAL*8 :: vx, vy, vz, idx2, &
              dxL_14, dxR_14, dxL_34, dxR_34, &
              dyL_14, dyR_14, dyL_34, dyR_34, &
              dzL_14, dzR_14, dzL_34, dzR_34, &
              dx_14_12, dx_34_12, dy_14_12, dy_34_12, dy_12_14, dy_12_34, dz_12_14, dz_12_34            

               
    idx2 = 2.d0 / dx
 
    ! Corner electric field z-component: Ez_i+1/2,j+1/2,k    
    DO k = klow-nb, khi+nb-1
        DO j = jlow-nb, jhi+nb-1
            DO i = ilow-nb, ihi+nb-1  
        
                ! get normal velocity at x-interface (roe-average value)
                vx = 0.5d0 * (q_3d(i,j,k,2)/q_3d(i,j,k,1) + q_3d(i+1,j,k,2)/q_3d(i+1,j,k,1) )
            
                ! get normal velocity at y-interface (roe-average value)
                vy = 0.5d0 * (q_3d(i,j,k,3)/q_3d(i,j,k,1) + q_3d(i,j+1,k,3)/q_3d(i,j+1,k,1) )
    
                ! compute upwinded interface gradients of electric field at x-interface
                dyL_14 = 2.d0 * (emf_3d_in(i,j,k,3,2) - emf_3d_in(i,j,k,3,3)) 
                dyR_14 = 2.d0 * (emf_3d_in(i+1,j,k,3,2) - emf_3d_in(i+1,j,k,3,3)) 
                dyL_34 = 2.d0 * (emf_3d_in(i,j+1,k,3,3) - emf_3d_in(i,j,k,3,2)) 
                dyR_34 = 2.d0 * (emf_3d_in(i+1,j+1,k,3,3) - emf_3d_in(i+1,j,k,3,2))

                IF(vx .GT. ZERO) THEN         
                    dy_12_14 =  dyL_14
                    dy_12_34 =  dyL_34                
                ELSE IF(vx .LT. -ZERO) THEN
                    dy_12_14 =  dyR_14
                    dy_12_34 =  dyR_34 
                ELSE
                    dy_12_14 = 0.5d0 * (dyL_14 + dyR_14)
                    dy_12_34 = 0.5d0 * (dyL_34 + dyR_34)
                END IF
            
                ! compute upwinded interface gradients of electric field at y-interface
                dxL_14 = 2.d0 * (emf_3d_in(i,j,k,3,1) - emf_3d_in(i,j,k,3,3))             
                dxR_14 = 2.d0 * (emf_3d_in(i,j+1,k,3,1) - emf_3d_in(i,j+1,k,3,3))            
                dxL_34 = 2.d0 * (emf_3d_in(i+1,j,k,3,3) - emf_3d_in(i,j,k,3,1))
                dxR_34 = 2.d0 * (emf_3d_in(i+1,j+1,k,3,3) - emf_3d_in(i,j+1,k,3,1))
                        
                IF(vy .GT. ZERO) THEN         
                    dx_14_12 =  dxL_14
                    dx_34_12 =  dxL_34 
                ELSE IF(vy .LT. -ZERO) THEN
                    dx_14_12 =  dxR_14
                    dx_34_12 =  dxR_34 
                ELSE
                    dx_14_12 = 0.5d0 * (dxL_14 + dxR_14)
                    dx_34_12 = 0.5d0 * (dxL_34 + dxR_34)
                END IF
             
                ! upwinded cell-corner electric field z component: Ez_i+1/2,j+1/2,k 
                emf_corner_in(i,j,k,3) =  0.25d0 * (emf_3d_in(i,j,k,3,1) + emf_3d_in(i,j+1,k,3,1) + &
                                       emf_3d_in(i,j,k,3,2) + emf_3d_in(i+1,j,k,3,2)) + &
                                       0.125d0 * ((dy_12_14 - dy_12_34) + (dx_14_12 - dx_34_12)) 
               
                !######## non-upwinded version##########
                !emf_corner_in(i,j,k,3) =  0.25d0 * (emf_3d_in(i,j,k,3,1) + emf_3d_in(i,j+1,k,3,1) + &
                !                       emf_3d_in(i,j,k,3,2) + emf_3d_in(i+1,j,k,3,2))                

                 
            END DO
        END DO    
    END DO    


    ! Corner electric field y-component: Ey_i+1/2,j,k+1/2
    DO k = klow-nb, khi+nb-1
        DO j = jlow-nb, jhi+nb-1
            DO i = ilow-nb, ihi+nb-1  
        
                ! get normal velocity at x-interface (roe-average value)
                vx = 0.5d0 * (q_3d(i,j,k,2)/q_3d(i,j,k,1) + q_3d(i+1,j,k,2)/q_3d(i+1,j,k,1) )
            
                ! get normal velocity at z-interface (roe-average value)
                vz = 0.5d0 * (q_3d(i,j,k,4)/q_3d(i,j,k,1) + q_3d(i,j,k+1,4)/q_3d(i,j,k+1,1) )
    
                ! compute upwinded interface gradients of electric field at x-interface
                dzL_14 = 2.d0 * (emf_3d_in(i,j,k,2,2) - emf_3d_in(i,j,k,2,3)) 
                dzR_14 = 2.d0 * (emf_3d_in(i+1,j,k,2,2) - emf_3d_in(i+1,j,k,2,3)) 
                dzL_34 = 2.d0 * (emf_3d_in(i,j,k+1,2,3) - emf_3d_in(i,j,k,2,2)) 
                dzR_34 = 2.d0 * (emf_3d_in(i+1,j,k+1,2,3) - emf_3d_in(i+1,j,k,2,2))
                
                IF(vx .GT. ZERO) THEN         
                    dz_12_14 =  dzL_14
                    dz_12_34 =  dzL_34                
                ELSE IF(vx .LT. -ZERO) THEN
                    dz_12_14 =  dzR_14
                    dz_12_34 =  dzR_34 
                ELSE
                    dz_12_14 = 0.5d0 * (dzL_14 + dzR_14)
                    dz_12_34 = 0.5d0 * (dzL_34 + dzR_34)
                END IF
            
                ! compute upwinded interface gradients of electric field at z-interface
                dxL_14 = 2.d0 * (emf_3d_in(i,j,k,2,1) - emf_3d_in(i,j,k,2,3))             
                dxR_14 = 2.d0 * (emf_3d_in(i,j,k+1,2,1) - emf_3d_in(i,j,k+1,2,3))            
                dxL_34 = 2.d0 * (emf_3d_in(i+1,j,k,2,3) - emf_3d_in(i,j,k,2,1))
                dxR_34 = 2.d0 * (emf_3d_in(i+1,j,k+1,2,3) - emf_3d_in(i,j,k+1,2,1))
                        
                IF(vz .GT. ZERO) THEN         
                    dx_14_12 =  dxL_14
                    dx_34_12 =  dxL_34 
                ELSE IF(vz .LT. -ZERO) THEN
                    dx_14_12 =  dxR_14
                    dx_34_12 =  dxR_34 
                ELSE
                    dx_14_12 = 0.5d0 * (dxL_14 + dxR_14)
                    dx_34_12 = 0.5d0 * (dxL_34 + dxR_34)
                END IF
             
                ! upwinded cell-corner electric field y component: Ey_i+1/2,j,k+1/2 
                emf_corner_in(i,j,k,2) =  0.25d0 * (emf_3d_in(i,j,k,2,1) + emf_3d_in(i,j,k+1,2,1) + &
                                       emf_3d_in(i,j,k,2,2) + emf_3d_in(i+1,j,k,2,2)) + &
                                       0.125d0 * ((dz_12_14 - dz_12_34) + (dx_14_12 - dx_34_12)) 
                 
                !######## non-upwinded version##########
                !emf_corner_in(i,j,k,2) =  0.25d0 * (emf_3d_in(i,j,k,2,1) + emf_3d_in(i,j,k+1,2,1) + &
                !                       emf_3d_in(i,j,k,2,2) + emf_3d_in(i+1,j,k,2,2)) 
                                       
            END DO
        END DO    
    END DO    

    ! Corner electric field x-component: Ex_i,j+1/2,k+1/2
    DO k = klow-nb, khi+nb-1
        DO j = jlow-nb, jhi+nb-1
            DO i = ilow-nb, ihi+nb-1  
        
                ! get normal velocity at y-interface (roe-average value)
                vy = 0.5d0 * (q_3d(i,j,k,3)/q_3d(i,j,k,1) + q_3d(i,j+1,k,3)/q_3d(i,j+1,k,1) )
                
                ! get normal velocity at z-interface (roe-average value)
                vz = 0.5d0 * (q_3d(i,j,k,4)/q_3d(i,j,k,1) + q_3d(i,j,k+1,4)/q_3d(i,j,k+1,1) )
    
                ! compute upwinded interface gradients of electric field at y-interface
                dzL_14 = 2.d0 * (emf_3d_in(i,j,k,1,2) - emf_3d_in(i,j,k,1,3)) 
                dzR_14 = 2.d0 * (emf_3d_in(i,j+1,k,1,2) - emf_3d_in(i,j+1,k,1,3)) 
                dzL_34 = 2.d0 * (emf_3d_in(i,j,k+1,1,3) - emf_3d_in(i,j,k,1,2)) 
                dzR_34 = 2.d0 * (emf_3d_in(i,j+1,k+1,1,3) - emf_3d_in(i,j+1,k,1,2))
                
                IF(vy .GT. ZERO) THEN         
                    dz_12_14 =  dzL_14
                    dz_12_34 =  dzL_34                
                ELSE IF(vy .LT. -ZERO) THEN
                    dz_12_14 =  dzR_14
                    dz_12_34 =  dzR_34 
                ELSE
                    dz_12_14 = 0.5d0 * (dzL_14 + dzR_14)
                    dz_12_34 = 0.5d0 * (dzL_34 + dzR_34)
                END IF
            
                ! compute upwinded interface gradients of electric field at z-interface
                dyL_14 = 2.d0 * (emf_3d_in(i,j,k,1,1) - emf_3d_in(i,j,k,1,3))             
                dyR_14 = 2.d0 * (emf_3d_in(i,j,k+1,1,1) - emf_3d_in(i,j,k+1,1,3))            
                dyL_34 = 2.d0 * (emf_3d_in(i,j+1,k,1,3) - emf_3d_in(i,j,k,1,1))
                dyR_34 = 2.d0 * (emf_3d_in(i,j+1,k+1,1,3) - emf_3d_in(i,j,k+1,1,1))
                        
                IF(vz .GT. ZERO) THEN         
                    dy_14_12 =  dyL_14
                    dy_34_12 =  dyL_34 
                ELSE IF(vz .LT. -ZERO) THEN
                    dy_14_12 =  dyR_14
                    dy_34_12 =  dyR_34 
                ELSE
                    dy_14_12 = 0.5d0 * (dyL_14 + dyR_14)
                    dy_34_12 = 0.5d0 * (dyL_34 + dyR_34)
                END IF
             
                ! upwinded cell-corner electric field x component: Ex_i,j+1/2,k+1/2 
                emf_corner_in(i,j,k,1) =  0.25d0 * (emf_3d_in(i,j,k,1,1) + emf_3d_in(i,j,k+1,1,1) + &
                                       emf_3d_in(i,j,k,1,2) + emf_3d_in(i,j+1,k,1,2)) + &
                                       0.125d0 * ((dz_12_14 - dz_12_34) + (dy_14_12 - dy_34_12)) 
                
                !######## non-upwinded version##########                
                !emf_corner_in(i,j,k,1) =  0.25d0 * (emf_3d_in(i,j,k,1,1) + emf_3d_in(i,j,k+1,1,1) + &
                !                       emf_3d_in(i,j,k,1,2) + emf_3d_in(i,j+1,k,1,2)) 
                 
            END DO
        END DO    
    END DO    


END SUBROUTINE compute_cell_corner_emf_3d


! final (unsplit) finite volume update of the state vector (excluding Bx and By, which require CT update)
SUBROUTINE update_hydro_state_3d(dtdx, ilow, ihi, jlow, jhi, klow, khi)

    REAL*8, INTENT(IN) :: dtdx
    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi, klow, khi
    INTEGER :: i, j, k
  
    IF(flux_protection) THEN

        protection_flags = .FALSE.
   
        ! make a copy of the initial state
     DO k = klow-nb, khi+nb    
         DO j = jlow-nb, jhi+nb     
             DO i = ilow-nb, ihi+nb
                    q_3d_prot(i,j,k,:) = q_3d(i,j,k,:)
                END DO
            END DO
        END DO     
    
    END IF
    
    ! update the conserved variables (density and momentum only)
     DO k = klow-nb+1, khi+nb-1    
         DO j = jlow-nb+1, jhi+nb-1     
             DO i = ilow-nb+1, ihi+nb-1
        
                q_3d(i,j,k,1) = q_3d(i,j,k,1) - dtdx * ( (Fs(i,j,k,1,1) - Fs(i-1,j,k,1,1)) + (Fs(i,j,k,1,2) - Fs(i,j-1,k,1,2)) + &
                                (Fs(i,j,k,1,3) - Fs(i,j,k-1,1,3)) ) 
                                
                q_3d(i,j,k,2) = q_3d(i,j,k,2) - dtdx * ( (Fs(i,j,k,2,1) - Fs(i-1,j,k,2,1)) + (Fs(i,j,k,2,2) - Fs(i,j-1,k,2,2)) + &
                                (Fs(i,j,k,2,3) - Fs(i,j,k-1,2,3)) )
                                
                q_3d(i,j,k,3) = q_3d(i,j,k,3) - dtdx * ( (Fs(i,j,k,3,1) - Fs(i-1,j,k,3,1)) + (Fs(i,j,k,3,2) - Fs(i,j-1,k,3,2)) + &
                                (Fs(i,j,k,3,3) - Fs(i,j,k-1,3,3)) )
                                
                q_3d(i,j,k,4) = q_3d(i,j,k,4) - dtdx * ( (Fs(i,j,k,4,1) - Fs(i-1,j,k,4,1)) + (Fs(i,j,k,4,2) - Fs(i,j-1,k,4,2)) + &
                                (Fs(i,j,k,4,3) - Fs(i,j,k-1,4,3)) )
                       
                !IF(CT_off) THEN
                !
                !    q_3d(i,j,k,5) = q_3d(i,j,k,5) - dtdx * ( (Fs(i,j,k,5,1) - Fs(i-1,j,k,5,1)) + (Fs(i,j,k,5,2) - Fs(i,j-1,k,5,2)) + &
                !                    (Fs(i,j,k,5,3) - Fs(i,j,k-1,5,3)) )  
                !                        
                !    q_3d(i,j,k,6) = q_3d(i,j,k,6) - dtdx * ( (Fs(i,j,k,6,1) - Fs(i-1,j,k,6,1)) + (Fs(i,j,k,6,2) - Fs(i,j-1,k,6,2)) + &
                !                    (Fs(i,j,k,6,3) - Fs(i,j,k-1,6,3)) )  
                !                        
                !    q_3d(i,j,k,7) = q_3d(i,j,k,7) - dtdx * ( (Fs(i,j,k,7,1) - Fs(i-1,j,k,7,1)) + (Fs(i,j,k,7,2) - Fs(i,j-1,k,7,2)) + &
                !                    (Fs(i,j,k,7,3) - Fs(i,j,k-1,7,3)) )
                !
                !END IF       
                       
            END DO    
        END DO
    END DO

    ! check for unphysical state
    IF(flux_protection) THEN

     DO k = klow-nb+1, khi+nb-1    
         DO j = jlow-nb+1, jhi+nb-1     
             DO i = ilow-nb+1, ihi+nb-1
       
                    IF(q_3d(i,j,k,1) .LT. MIN_DENS) THEN
                    
                        protection_flags(i,j,k) = .TRUE.
                        
                        ! copy protection fluxes through boundary of bad cell into main Fs array
                        Fs(i,j,k,:,1)   = Fs_prot(i,j,k,:,1)
                        Fs(i-1,j,k,:,1) = Fs_prot(i-1,j,k,:,1)
                        Fs(i,j,k,:,2)   = Fs_prot(i,j,k,:,2)
                        Fs(i,j-1,k,:,2) = Fs_prot(i,j-1,k,:,2)
                        Fs(i,j,k,:,3)   = Fs_prot(i,j,k,:,3)
                        Fs(i,j,k-1,:,3) = Fs_prot(i,j,k-1,:,3)
                        
                        ! also copy the protection cell corner emf so that they will automatically get
                        ! applied during the CT update
                        emf_corner(i,j,k,1)     =  emf_corner_prot(i,j,k,1)
                        emf_corner(i-1,j,k,1)   =  emf_corner_prot(i-1,j,k,1)
                        emf_corner(i,j-1,k,1)   =  emf_corner_prot(i,j-1,k,1)
                        emf_corner(i,j,k-1,1)   =  emf_corner_prot(i,j,k-1,1)
                        
                        emf_corner(i,j,k,2)     =  emf_corner_prot(i,j,k,2)
                        emf_corner(i-1,j,k,2)   =  emf_corner_prot(i-1,j,k,2)
                        emf_corner(i,j-1,k,2)   =  emf_corner_prot(i,j-1,k,2)
                        emf_corner(i,j,k-1,2)   =  emf_corner_prot(i,j,k-1,2)
                        
                        emf_corner(i,j,k,3)     =  emf_corner_prot(i,j,k,3)
                        emf_corner(i-1,j,k,3)   =  emf_corner_prot(i-1,j,k,3)
                        emf_corner(i,j-1,k,3)   =  emf_corner_prot(i,j-1,k,3)
                        emf_corner(i,j,k-1,3)   =  emf_corner_prot(i,j,k-1,3)
                        
                    END IF
       
                END DO
            END DO
        END DO

     DO k = klow-nb+1, khi+nb-1    
         DO j = jlow-nb+1, jhi+nb-1     
             DO i = ilow-nb+1, ihi+nb-1

                     IF(protection_flags(i,j,k) .OR. protection_flags(i-1,j,k) .OR. protection_flags(i+1,j,k) .OR. &
                        protection_flags(i,j-1,k) .OR. protection_flags(i,j+1,k) .OR. protection_flags(i,j,k-1) .OR. protection_flags(i,j,k+1)) THEN      

                        q_3d(i,j,k,1) = q_3d_prot(i,j,k,1) - dtdx * ( (Fs(i,j,k,1,1) - Fs(i-1,j,k,1,1)) + (Fs(i,j,k,1,2) - Fs(i,j-1,k,1,2)) + &
                                        (Fs(i,j,k,1,3) - Fs(i,j,k-1,1,3)) ) 
                                        
                        q_3d(i,j,k,2) = q_3d_prot(i,j,k,2) - dtdx * ( (Fs(i,j,k,2,1) - Fs(i-1,j,k,2,1)) + (Fs(i,j,k,2,2) - Fs(i,j-1,k,2,2)) + &
                                        (Fs(i,j,k,2,3) - Fs(i,j,k-1,2,3)) )
                                        
                        q_3d(i,j,k,3) = q_3d_prot(i,j,k,3) - dtdx * ( (Fs(i,j,k,3,1) - Fs(i-1,j,k,3,1)) + (Fs(i,j,k,3,2) - Fs(i,j-1,k,3,2)) + &
                                        (Fs(i,j,k,3,3) - Fs(i,j,k-1,3,3)) )
                                        
                        q_3d(i,j,k,4) = q_3d_prot(i,j,k,4) - dtdx * ( (Fs(i,j,k,4,1) - Fs(i-1,j,k,4,1)) + (Fs(i,j,k,4,2) - Fs(i,j-1,k,4,2)) + &
                                        (Fs(i,j,k,4,3) - Fs(i,j,k-1,4,3)) )
                                      
                    END IF
       
                END DO
            END DO
        END DO

    END IF

    CALL final_check_unphysical(ilow-nb+1, ihi+nb-1, jlow-nb+1, jhi+nb-1, klow-nb+1, khi+nb-1)


END SUBROUTINE update_hydro_state_3d


! updates cell-faces magnetic fields via CT
SUBROUTINE update_CT_3d(dtdx, ilow, ihi, jlow, jhi, klow, khi)

    REAL*8, INTENT(IN) :: dtdx
    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi, klow, khi
    INTEGER :: i, j, k
    
    
    DO k = klow-nb+1, khi+nb-1
        DO j = jlow-nb+1, jhi+nb-1
            DO i = ilow-nb+1, ihi+nb-1  

                ! Bx_i+1/2,j,k  
                bface_3d(i,j,k,1) = bface_3d(i,j,k,1) - dtdx * ( (emf_corner(i,j,k,3) - emf_corner(i,j-1,k,3)) - &
                                    (emf_corner(i,j,k,2) - emf_corner(i,j,k-1,2)) )               
               
                ! By_i,j+1/2,k
                bface_3d(i,j,k,2) = bface_3d(i,j,k,2) - dtdx * ( (emf_corner(i,j,k,1) - emf_corner(i,j,k-1,1)) - &
                                    (emf_corner(i,j,k,3) - emf_corner(i-1,j,k,3)) )          

                ! Bz_i,j,k+1/2
                bface_3d(i,j,k,3) = bface_3d(i,j,k,3) - dtdx * ( (emf_corner(i,j,k,2) - emf_corner(i-1,j,k,2)) - &
                                    (emf_corner(i,j,k,1) - emf_corner(i,j-1,k,1)) )                                              
                                                   
            END DO
        END DO
    END DO


    ! Also compute the cell-centered magnetic fields
    ! by averaging over the cell-face values (i.e. linear interpolation)    
    DO k = klow-nb+2, khi+nb-1
        DO j = jlow-nb+2, jhi+nb-1
            DO i = ilow-nb+2, ihi+nb-1  
  
                ! Bx_i,j,k
                q_3d(i,j,k,5) = 0.5d0 * (bface_3d(i,j,k,1) + bface_3d(i-1,j,k,1))  
                
                ! By_i,j,k
                q_3d(i,j,k,6) = 0.5d0 * (bface_3d(i,j,k,2) + bface_3d(i,j-1,k,2))

                ! Bz_i,j,k
                q_3d(i,j,k,7) = 0.5d0 * (bface_3d(i,j,k,3) + bface_3d(i,j,k-1,3))                    
               
            END DO
        END DO
    END DO
 
    
END SUBROUTINE update_CT_3d



! Compute cell-interface roe averaged state
SUBROUTINE compute_roe_avg_state_1d(ilow, ihi, jlow, jhi)

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi
    REAL*8 :: rho2, vx2, vy2, vz2, Bx2, By2, Bz2, isrho
    INTEGER :: i, j
    
       
    DO j = jlow-nb, jhi+nb     
        DO i = ilow-nb, ihi+nb-1

            ! compute Roe-averaged primitive variables (this is just an arithmetic average)
            rho2 = 0.5d0 * ( qintL(i,j,1) + qintR(i,j,1) )
            isrho = 1.D0 / DSQRT(rho2)
            vx2  = 0.5d0 * ( qintL(i,j,2)/qintL(i,j,1) + qintR(i,j,2)/qintR(i,j,1) )
            vy2  = 0.5d0 * ( qintL(i,j,3)/qintL(i,j,1) + qintR(i,j,3)/qintR(i,j,1) )
            vz2  = 0.5d0 * ( qintL(i,j,4)/qintL(i,j,1) + qintR(i,j,4)/qintR(i,j,1) )
            Bx2  = qintL(i,j,5)  ! OR qintR(i,j,5), same thing
            By2  = 0.5d0 * ( qintL(i,j,6) + qintR(i,j,6) )
            Bz2  = 0.5d0 * ( qintL(i,j,7) + qintR(i,j,7) )
            
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
SUBROUTINE compute_roe_eigenvalues_1d(ilow, ihi, jlow, jhi)

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi
    REAL*8 :: rho2, vx2, vy2, vz2, bx2, by2, bz2, asqr
    REAL*8 :: cf, ca, cs, bsqr, btsqr, bx, tmp1, tmp2, isrho
    INTEGER :: i, j
        

    asqr  = sound_speed**2


    DO j = jlow-nb, jhi+nb     
        DO i = ilow-nb, ihi+nb-1

            rho2  = qavg(i,j,1)
            vx2   = qavg(i,j,2)
            bx    = qavg(i,j,5)
            bsqr  = qavg(i,j,5)**2 + qavg(i,j,6)**2 + qavg(i,j,7)**2
            btsqr =  qavg(i,j,6)**2 + qavg(i,j,7)**2
            tmp1  = asqr + bsqr        
            tmp2  = DSQRT( (asqr - bsqr)**2 + 4.d0 * asqr * btsqr )
            !tmp2  = DSQRT(MAX(0.D0, tmp1**2 - 4.d0 * asqr * bx**2) )

            cs = DSQRT( MAX(0.D0, 0.5d0 * (tmp1 - tmp2) ))  ! slow mode speed
            ca = ABS(bx)                                    ! alfven speed
            cf = DSQRT(0.5d0 * (tmp1 + tmp2))               ! fast mode speed
            
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

            ! check max wave speed
            smax = MAX(smax, ABS(vx2)+cf)


            !IF(cf .LT. ca) THEN
            !    PRINT*,'Round off error. cf - ca = ',cf-ca
            !    STOP        
            !END IF

        END DO 
    END DO
    
    
END SUBROUTINE compute_roe_eigenvalues_1d


! Compute eigenvectors
SUBROUTINE compute_roe_eigenvectors_1d(ilow, ihi, jlow, jhi)

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi
    REAL*8 :: rho2, vx2, vy2, vz2, bx2, by2, bz2, a1, asqr, isrho, amc
    REAL*8 :: cf, ca, cs, bsqr, bt, tmp1, tmp2, tmp3, tmp4, tmp5, signbx, signbt, itheta1, itheta2
    REAL*8 :: alphas(1-nb:nmax+nb,1-nb:nmax+nb), alphaf(1-nb:nmax+nb,1-nb:nmax+nb), &
              betay(1-nb:nmax+nb,1-nb:nmax+nb), betaz(1-nb:nmax+nb,1-nb:nmax+nb)
    REAL*8 :: delq(6), Lminus(6), Lplus(6)
    INTEGER :: i, j, k 
     

    asqr   = sound_speed**2      

         
    DO j = jlow-nb, jhi+nb     
        DO i = ilow-nb, ihi+nb-1

            bt = DSQRT(qavg(i,j,6)**2 + qavg(i,j,7)**2) ! magnitude of B-transverse

            ! compute eigenvector renormalization co-efficients, taking limiting values for degenracies (K98 eqn. 2.13-15)
            
            !*************************** ORIGINAL *************************************
            !  Fails if c_s = c_a = c_f = c_sound
            !
            !IF(ABS(wave_speeds(i,j,3)**2 - wave_speeds(i,j,1)**2) .LT. ZERO) THEN
            !    alphaf(i,j) = 1.D0
            !    alphas(i,j) = 1.D0
            !ELSE
            !    tmp1 = 1.D0 / DSQRT( MAX(0.d0, wave_speeds(i,j,3)**2 - wave_speeds(i,j,1)**2) )
            !    alphaf(i,j) = DSQRT( MAX(0.d0, wave_speeds(i,j,3)**2 - wave_speeds(i,j,2)**2) ) * tmp1      ! Need to be careful with square roots. Sometimes round-off errors can make cf < ca, or ca < cs.   
            !    alphas(i,j) = DSQRT( MAX(0.d0, wave_speeds(i,j,3)**2 - asqr) ) * tmp1
            !END IF
            
            
            !************************** MODIFIED V1 ********************************        
            IF(ABS(wave_speeds(i,j,3)**2 - wave_speeds(i,j,1)**2) .LT. ZERO) THEN
                alphaf(i,j) = 1.D0
                alphas(i,j) = 1.D0
            ELSE IF(ABS(wave_speeds(i,j,3)**2 - wave_speeds(i,j,2)**2) .LT. ZERO) THEN
                alphaf(i,j) = 0.D0
                alphas(i,j) = 1.D0
            ELSE IF(ABS(wave_speeds(i,j,3)**2 - asqr) .LT. ZERO) THEN
                alphaf(i,j) = 1.D0
                alphas(i,j) = 0.D0            
            ELSE
                tmp1 = 1.D0 / DSQRT( MAX(0.d0, wave_speeds(i,j,3)**2 - wave_speeds(i,j,1)**2) )
                alphaf(i,j) = DSQRT( MAX(0.d0, wave_speeds(i,j,3)**2 - wave_speeds(i,j,2)**2) ) * tmp1      ! Need to be careful with square roots. Sometimes round-off errors can make cf < ca, or ca < cs.   
                alphas(i,j) = DSQRT( MAX(0.d0, wave_speeds(i,j,3)**2 - asqr) ) * tmp1
            END IF
            !*************************************************************************
             
         
            !************************** MODIFIED V2 ********************************        
            !IF((bt .LT. ZERO .AND. ABS(wave_speeds(i,j,2)**2 - asqr) .LT. ZERO) .OR. &
            !   (ABS(wave_speeds(i,j,3)**2 - wave_speeds(i,j,1)**2) .LT. ZERO) ) THEN        
            !    alphaf(i,j) = 1.D0
            !    alphas(i,j) = 1.D0   
            !ELSE
            !    tmp1 = 1.D0 / DSQRT( MAX(0.d0, wave_speeds(i,j,3)**2 - wave_speeds(i,j,1)**2) )
            !    alphaf(i,j) = DSQRT( MAX(0.d0, wave_speeds(i,j,3)**2 - wave_speeds(i,j,2)**2) ) * tmp1      ! Need to be careful with square roots. Sometimes round-off errors can make cf < ca, or ca < cs.   
            !    alphas(i,j) = DSQRT( MAX(0.d0, wave_speeds(i,j,3)**2 - asqr) ) * tmp1 
            !END IF
            !*************************************************************************


            
            !IF(alphaf(i,j) .LT. ZERO .AND. alphas(i,j) .LT. ZERO) THEN
            !    PRINT*,'alphaf, alphas =',alphaf(i,j), alphas(i,j)
            !    PRINT*,'bt, ca^2-a^2, cf^2-cs^2 = ', bt, wave_speeds(i,j,2)**2 - asqr, wave_speeds(i,j,3)**2 - wave_speeds(i,j,1)**2
            !    STOP
            !END IF


            IF(bt .LT. ZERO) THEN        
                betay(i,j) = 1.D0 / DSQRT(2.D0)
                betaz(i,j) = 1.D0 / DSQRT(2.D0)  
            ELSE
                tmp2 = 1.D0 / bt
                betay(i,j) = qavg(i,j,6) * tmp2
                betaz(i,j) = qavg(i,j,7) * tmp2       
            END IF
            
        END DO    
    END DO
    
    !*************************************************
    ! compute the right eigenvectors (K98 eq 2.12)
    !*************************************************
    
    DO j = jlow-nb, jhi+nb     
        DO i = ilow-nb, ihi+nb-1
        
            isrho = 1.D0 / DSQRT(qavg(i,j,1))       
            bt = DSQRT(qavg(i,j,6)**2 + qavg(i,j,7)**2) ! magnitude of B-transverse
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

    DO j = jlow-nb, jhi+nb     
        DO i = ilow-nb, ihi+nb-1
        
            isrho = 1.D0 / DSQRT(qavg(i,j,1))       
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
    
    DO j = jlow-nb, jhi+nb     
        DO i = ilow-nb, ihi+nb-1
        
            isrho = 1.D0 / DSQRT(qavg(i,j,1))       
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
        
    !*****************************************************
    ! compute the characteristic variables (K98 eqn. 2.27)
    !*****************************************************
    
    DO j = jlow-nb, jhi+nb     
        DO i = ilow-nb, ihi+nb-1
                                
            signbx = SIGN(1.D0, qavg(i,j,5))
            amc = sound_speed**2 - wave_speeds(i,j,2)**2

            ! q_i+1 - q_i (Note: The Bx component has to be excluded)
            delq(1:4) = qintR(i,j,1:4) - qintL(i,j,1:4) 
            delq(5)   = qintR(i,j,6)   - qintL(i,j,6)
            delq(6)   = qintR(i,j,7)   - qintL(i,j,7)
        
            itheta1 = 1.D0 / (2.D0 * ( alphaf(i,j)**2 * sound_speed**2 + alphas(i,j)**2 * wave_speeds(i,j,3)**2))  
            itheta2 = 1.D0 / (2.D0 * ( alphaf(i,j)**2 * sound_speed * wave_speeds(i,j,3) + &
                      alphas(i,j)**2 * wave_speeds(i,j,1) * wave_speeds(i,j,2) ) )
           
            tmp1 = itheta1 * alphaf(i,j) * sound_speed**2
            tmp2 = itheta2 * ( - alphaf(i,j) * sound_speed * qavg(i,j,2) + &
                   alphas(i,j) * wave_speeds(i,j,1) * signbx * ( betay(i,j) * qavg(i,j,3) + betaz(i,j) * qavg(i,j,4)) )         
     
            tmp3 = itheta2 * alphas(i,j) * wave_speeds(i,j,1) * signbx 
            tmp4 = itheta1 * alphas(i,j) * wave_speeds(i,j,3) * DSQRT(qavg(i,j,1))    
            
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
    

    DO j = jlow-nb, jhi+nb     
        DO i = ilow-nb, ihi+nb-1
                  
            signbx = SIGN(1.D0, qavg(i,j,5))

            ! q_i+1 - q_i (Note: The Bx component has to be excluded)
            delq(1:4) = qintR(i,j,1:4) - qintL(i,j,1:4) 
            delq(5)   = qintR(i,j,6)   - qintL(i,j,6)
            delq(6)   = qintR(i,j,7)   - qintL(i,j,7)
            
            tmp1 = 0.5d0 * (betaz(i,j) * qavg(i,j,3) - betay(i,j) * qavg(i,j,4)) * signbx
            tmp2 = 0.5d0 * signbx
            tmp3 = 0.5d0 * DSQRT(qavg(i,j,1))
            
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
    
        
    DO j = jlow-nb, jhi+nb     
        DO i = ilow-nb, ihi+nb-1
                  
            signbx = SIGN(1.D0, qavg(i,j,5))
            amc = sound_speed**2 - wave_speeds(i,j,2)**2 

            ! q_i+1 - q_i (Note: The Bx component has to be excluded)
            delq(1:4) = qintR(i,j,1:4) - qintL(i,j,1:4) 
            delq(5)   = qintR(i,j,6)   - qintL(i,j,6)
            delq(6)   = qintR(i,j,7)   - qintL(i,j,7)
            
            itheta1 = 1.D0 / (2.D0 * ( alphaf(i,j)**2 * sound_speed**2 + alphas(i,j)**2 * wave_speeds(i,j,3)**2))  
            itheta2 = 1.D0 / (2.D0 * ( alphaf(i,j)**2 * sound_speed * wave_speeds(i,j,3) + &
                      alphas(i,j)**2 * wave_speeds(i,j,1) * wave_speeds(i,j,2) ) )
           
            tmp1 = itheta1 * alphas(i,j) * wave_speeds(i,j,3)**2
            tmp2 = itheta2 * ( alphas(i,j) * wave_speeds(i,j,2) * qavg(i,j,2) + &
                   alphaf(i,j) * wave_speeds(i,j,3) * signbx * ( betay(i,j) * qavg(i,j,3) + betaz(i,j) * qavg(i,j,4)) )         
     
            tmp3 = itheta2 * alphaf(i,j) * wave_speeds(i,j,3) * signbx 
            tmp4 = itheta1 * alphaf(i,j) * wave_speeds(i,j,3) * DSQRT(qavg(i,j,1))    
            
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
        

END SUBROUTINE compute_roe_eigenvectors_1d


! Compute the time-averaged cell-interface fluxes using Roe's riemann solver
SUBROUTINE compute_ROE_fluxes_1d(Flux, ilow, ihi, jlow, jhi)

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi
    REAL*8, INTENT(INOUT) :: Flux(1-nb:nmax+nb,1-nb:nmax+nb,7)
    INTEGER :: i, j, k
    REAL*8 :: dxdt, gammak, betak
    REAL*8 :: F(1-nb:ihi+nb, 1-nb:jhi+nb,nwaves,2), FTVD(1-nb:ihi+nb, 1-nb:jhi+nb,nwaves)
    REAL*8 :: rho, vx, vy, vz, Bx, By, Bz

    
    ! clear temp array
    FTVD = 0.D0
    
    ! compute cell center fluxes 
    
    ! left fluxes
    DO j = jlow-nb, jhi+nb     
        DO i = ilow-nb, ihi+nb-1
        
            rho = qintL(i,j,1)
            vx  = qintL(i,j,2) / rho
            vy  = qintL(i,j,3) / rho
            vz  = qintL(i,j,4) / rho
            Bx  = qintL(i,j,5)
            By  = qintL(i,j,6)
            Bz  = qintL(i,j,7)
        
            F(i,j,1,1) = rho * vx
            F(i,j,2,1) = rho * (vx**2 + sound_speed**2) + 0.5d0 * (By**2 + Bz**2 - Bx**2)
            F(i,j,3,1) = rho * vx * vy - Bx * By
            F(i,j,4,1) = rho * vx * vz - Bx * Bz
            F(i,j,5,1) = vx * By - vy * Bx
            F(i,j,6,1) = vx * Bz - vz * Bx   
        
        END DO
    END DO

    ! right fluxes
    DO j = jlow-nb, jhi+nb     
        DO i = ilow-nb, ihi+nb-1
        
            rho = qintR(i,j,1)
            vx  = qintR(i,j,2) / rho
            vy  = qintR(i,j,3) / rho
            vz  = qintR(i,j,4) / rho
            Bx  = qintR(i,j,5)
            By  = qintR(i,j,6)
            Bz  = qintR(i,j,7)
        
            F(i,j,1,2) = rho * vx
            F(i,j,2,2) = rho * (vx**2 + sound_speed**2) + 0.5d0 * (By**2 + Bz**2 - Bx**2)
            F(i,j,3,2) = rho * vx * vy - Bx * By
            F(i,j,4,2) = rho * vx * vz - Bx * Bz
            F(i,j,5,2) = vx * By - vy * Bx
            F(i,j,6,2) = vx * Bz - vz * Bx   
        
        END DO
    END DO
    
   
    DO k = 1, nwaves
        DO j = jlow-nb, jhi+nb     
            DO i = ilow-nb, ihi+nb-1
            
            
                ! betak_i+1/2 
                
                !***************************************************
                ! without "entropy fix" (i.e. no extra dissipation)
                !***************************************************
                !betak = ABS(eigenvalues(i,j,k)) * ck(i,j,k)  
                
                !*****************************
                ! with Harten's "entropy fix" 
                !*****************************
                ! The goal of an entropy fix is to prevent the dissipation term (i.e. sum_k[betak * Rk(i,j,1,k)] )
                ! from approaching too close to zero. If the dissipation term becomes too small, it can cause smooth 
                ! parts of the solution to go bad.                 
                ! In Harten's scheme, a small extra amount to ABS(eigenvalue) if it becomes too small (i.e. smaller than some
                !  fixed amount delta). We can add different amounts to different wave modes depending on which modes need it more)
                IF(ABS(eigenvalues(i,j,k)) .GE. delta_h(k)) THEN
                    betak = ABS(eigenvalues(i,j,k)) * ck(i,j,k)  
                ELSE
                    betak = ((ABS(eigenvalues(i,j,k)) + delta_h(k)**2) / (2.d0 * delta_h(k))) * ck(i,j,k)  
                END IF
                
                ! accumulate flux contribution from each wave  
                FTVD(i,j,1) = FTVD(i,j,1) - betak * Rk(i,j,1,k) 
                FTVD(i,j,2) = FTVD(i,j,2) - betak * Rk(i,j,2,k) 
                FTVD(i,j,3) = FTVD(i,j,3) - betak * Rk(i,j,3,k) 
                FTVD(i,j,4) = FTVD(i,j,4) - betak * Rk(i,j,4,k) 
                FTVD(i,j,5) = FTVD(i,j,5) - betak * Rk(i,j,5,k) 
                FTVD(i,j,6) = FTVD(i,j,6) - betak * Rk(i,j,6,k) 
                
            END DO    
        END DO
    END DO

    ! add remaining term in the Roe flux
    DO j = jlow-nb, jhi+nb     
        DO i = ilow-nb, ihi+nb-1
            
            Flux(i,j,1) = 0.5d0 * ( FTVD(i,j,1) + F(i,j,1,1) + F(i,j,1,2) )
            Flux(i,j,2) = 0.5d0 * ( FTVD(i,j,2) + F(i,j,2,1) + F(i,j,2,2) )
            Flux(i,j,3) = 0.5d0 * ( FTVD(i,j,3) + F(i,j,3,1) + F(i,j,3,2) )
            Flux(i,j,4) = 0.5d0 * ( FTVD(i,j,4) + F(i,j,4,1) + F(i,j,4,2) )
            Flux(i,j,5) = 0.d0
            Flux(i,j,6) = 0.5d0 * ( FTVD(i,j,5) + F(i,j,5,1) + F(i,j,5,2) )
            Flux(i,j,7) = 0.5d0 * ( FTVD(i,j,6) + F(i,j,6,1) + F(i,j,6,2) )   
            
        END DO
    END DO
    
    
END SUBROUTINE compute_ROE_fluxes_1d


! This subroutine computes 1D "HLLE" fluxes
SUBROUTINE compute_HLLE_fluxes_1d(Flux, ilow, ihi, jlow, jhi)

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi
    REAL*8, INTENT(INOUT) :: Flux(1-nb:nmax+nb,1-nb:nmax+nb,7)
    INTEGER :: i, j, k
    REAL*8 :: uL(7), uR(7), Bx, vxL, vyL, vzL, vxR, vyR, vzR
    REAL*8 :: FL(7), FR(7)
    REAL*8 :: sL, sR, sL_s, sR_s, lambdaf_minus_L, lambdaf_minus_R, lambdaf_plus_L, lambdaf_plus_R
    REAL*8 :: tmp1, tmp2, asqr, bsqr, btsqr, cfL, cfR, isRL
    REAL*8 :: uavg(7), cf_avg, vx_avg, lambdaf_minus_avg, lambdaf_plus_avg
    REAL*8 :: uHLL(7), FHLL(7)    
    
    asqr = sound_speed**2
  
    DO j = jlow-nb, jhi+nb
        DO i = ilow-nb, ihi+nb-1
        
            ! left and right states and fluxes
            uL(:) = qintL(i,j,:)
            
            uR(:) = qintR(i,j,:)
            
            uavg(:) = 0.5d0 * (uL(:) + uR(:)) ! Roe-averaged state
            
            ! Need to make sure we're not feeding in bad interface states to the solver. Check for negative density           
            !uL(1) = MAX(MIN_DENS, uL(1))
            !uR(1) = MAX(MIN_DENS, uR(1))
            !uavg(1) = MAX(MIN_DENS, uavg(1))

            Bx = uL(5)  ! or uR(5), same thing..
          
            vxL = uL(2)/ uL(1)
            vyL = uL(3)/ uL(1)
            vzL = uL(4)/ uL(1)
            
            vxR = uR(2)/ uR(1)
            vyR = uR(3)/ uR(1)
            vzR = uR(4)/ uR(1)   
      
            FL(1) = uL(2)
            FL(2) = uL(1) * (vxL**2 + asqr) + 0.5d0 *(uL(6)**2 + uL(7)**2 - uL(5)**2)
            FL(3) = uL(2) * vyL - uL(5) * uL(6)
            FL(4) = uL(2) * vzL - uL(5) * uL(7)
            FL(5) = 0.d0
            FL(6) = vxL * uL(6) - vyL * uL(5)
            FL(7) = vxL * uL(7) - vzL * uL(5)   

            FR(1) = uR(2)
            FR(2) = uR(1) * (vxR**2 + asqr) + 0.5d0 *(uR(6)**2 + uR(7)**2 - uR(5)**2)
            FR(3) = uR(2) * vyR - uR(5) * uR(6)
            FR(4) = uR(2) * vzR - uR(5) * uR(7)
            FR(5) = 0.d0
            FR(6) = vxR * uR(6) - vyR * uR(5)
            FR(7) = vxR * uR(7) - vzR * uR(5)       
            
                        
            ! compute lower and upper signal speeds
            
            bsqr   = (uL(5)**2 + uL(6)**2 + uL(7)**2) / uL(1)
            btsqr  = (uL(6)**2 + uL(7)**2) / uL(1)
                        
            tmp1 = asqr + bsqr        
            tmp2 = DSQRT( (asqr - bsqr)**2 + 4.d0 * asqr * btsqr)
            cfL   = DSQRT( 0.5d0 * (tmp1 + tmp2) )
            
            lambdaf_minus_L = vxL - cfL  
            lambdaf_plus_L  = vxL + cfL  
            
            bsqr   = (uR(5)**2 + uR(6)**2 + uR(7)**2) / uR(1)
            btsqr  = (uR(6)**2 + uR(7)**2) / uR(1)
            
            tmp1 = asqr + bsqr        
            tmp2 = DSQRT( (asqr - bsqr)**2 + 4.d0 * asqr * btsqr)
            cfR   = DSQRT( 0.5d0 * (tmp1 + tmp2) )

            lambdaf_minus_R = vxR - cfR  
            lambdaf_plus_R  = vxR + cfR          


            bsqr  = (uavg(5)**2 + uavg(6)**2 + uavg(7)**2) / uavg(1)
            btsqr = (uavg(6)**2 + uavg(7)**2) / uavg(1)
            tmp1 = asqr + bsqr        
            tmp2 = DSQRT( (asqr - bsqr)**2 + 4.d0 * asqr * btsqr )
            cf_avg   = DSQRT( 0.5d0 * (tmp1 + tmp2) )
            
            vx_avg = uavg(2)/ uavg(1)

            lambdaf_minus_avg = vx_avg - cf_avg 
            lambdaf_plus_avg  = vx_avg + cf_avg

            sL = MIN(lambdaf_minus_L, lambdaf_minus_avg, 0.d0) !, 0.d0)
            sR = MAX(lambdaf_plus_R, lambdaf_plus_avg, 0.D0)
           
            ! check for max wave speed
            !smax = MAX(smax, ABS(sL), ABS(sR))
            smax = MAX(smax, ABS(vxL)+cfL, ABS(vxR)+cfR, ABS(vx_avg)+cf_avg)
            
            isRL = 1.d0 / (sR - sL)

            ! compute HLL state
            uHLL(:) = (sR * uR(:) - sL * uL(:) - FR(:) + FL(:) ) * isRL
             
            ! compute HLL flux
            FHLL(:)  = (sR * FL(:) - sL * FR(:) + sR * sL* (uR(:) - uL(:)) ) * isRL


            Flux(i,j,:) = FHLL(:)               
           
            
        END DO    
    END DO   
    
    
    
END SUBROUTINE compute_HLLE_fluxes_1d



! This subroutine computes 1D "HLLD" fluxes according to the algorithm from M07
SUBROUTINE compute_HLLD_fluxes_1d(Flux, ilow, ihi, jlow, jhi)

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi
    REAL*8, INTENT(INOUT) :: Flux(1-nb:nmax+nb,1-nb:nmax+nb,7)
    INTEGER :: i, j, k
    REAL*8 :: uL(7), uR(7), Bx, vxL, vyL, vzL, vxR, vyR, vzR
    REAL*8 :: FL(7), FR(7)
    REAL*8 :: rho_s, momx_s, vx_s, flux_rho_s, flux_momx_s, irho_s, isrho_s, X, iX 
    REAL*8 :: vy_s_L, vz_s_L, vy_s_R, vz_s_R, By_s_L, Bz_s_L , By_s_R, Bz_s_R
    REAL*8 :: vy_s_C, vz_s_C, By_s_C, Bz_s_C 
    REAL*8 :: sL, sR, sL_s, sR_s, lambdaf_minus_L, lambdaf_minus_R, lambdaf_plus_L, lambdaf_plus_R
    REAL*8 :: tmp1, tmp2, asqr, bsqr, btsqr, cfL, cfR, isRL
    REAL*8 :: uavg(7), cf_avg, vx_avg, lambdaf_minus_avg, lambdaf_plus_avg
    REAL*8 :: uHLL(7), FHLL(7)
    REAL*8 :: uL_s(7), uC_s(7), uR_s(7)
    
    asqr = sound_speed**2
  
   DO j = jlow-nb, jhi+nb
        DO i = ilow-nb, ihi+nb-1
             
        ! left and right states and fluxes
            uL(:) = qintL(i,j,:)
            
            uR(:) = qintR(i,j,:)
            
            uavg(:) = 0.5d0 * (uL(:) + uR(:)) ! Roe-averaged state
            
            ! Need to make sure we're not feeding in bad interface states to the solver. Check for negative density           
            !uL(1) = MAX(MIN_DENS, uL(1))
            !uR(1) = MAX(MIN_DENS, uR(1))
            !uavg(1) = MAX(MIN_DENS, uavg(1))

            Bx = uL(5)  ! or uR(5), same thing..
          
            vxL = uL(2)/ uL(1)
            vyL = uL(3)/ uL(1)
            vzL = uL(4)/ uL(1)
            
            vxR = uR(2)/ uR(1)
            vyR = uR(3)/ uR(1)
            vzR = uR(4)/ uR(1)   
      
            FL(1) = uL(2)
            FL(2) = uL(1) * (vxL**2 + asqr) + 0.5d0 *(uL(6)**2 + uL(7)**2 - uL(5)**2)
            FL(3) = uL(2) * vyL - uL(5) * uL(6)
            FL(4) = uL(2) * vzL - uL(5) * uL(7)
            FL(5) = 0.d0
            FL(6) = vxL * uL(6) - vyL * uL(5)
            FL(7) = vxL * uL(7) - vzL * uL(5)   

            FR(1) = uR(2)
            FR(2) = uR(1) * (vxR**2 + asqr) + 0.5d0 *(uR(6)**2 + uR(7)**2 - uR(5)**2)
            FR(3) = uR(2) * vyR - uR(5) * uR(6)
            FR(4) = uR(2) * vzR - uR(5) * uR(7)
            FR(5) = 0.d0
            FR(6) = vxR * uR(6) - vyR * uR(5)
            FR(7) = vxR * uR(7) - vzR * uR(5)       
            
                        
            ! compute lower and upper signal speeds
            
            bsqr   = (uL(5)**2 + uL(6)**2 + uL(7)**2) / uL(1)
            btsqr  = (uL(6)**2 + uL(7)**2) / uL(1)
                        
            tmp1 = asqr + bsqr        
            tmp2 = DSQRT( (asqr - bsqr)**2 + 4.d0 * asqr * btsqr)
            cfL   = DSQRT( 0.5d0 * (tmp1 + tmp2) )
            
            lambdaf_minus_L = vxL - cfL  
            lambdaf_plus_L  = vxL + cfL  
            
            bsqr   = (uR(5)**2 + uR(6)**2 + uR(7)**2) / uR(1)
            btsqr  = (uR(6)**2 + uR(7)**2) / uR(1)
            
            tmp1 = asqr + bsqr        
            tmp2 = DSQRT( (asqr - bsqr)**2 + 4.d0 * asqr * btsqr)
            cfR   = DSQRT( 0.5d0 * (tmp1 + tmp2) )

            lambdaf_minus_R = vxR - cfR  
            lambdaf_plus_R  = vxR + cfR          


            bsqr  = (uavg(5)**2 + uavg(6)**2 + uavg(7)**2) / uavg(1)
            btsqr = (uavg(6)**2 + uavg(7)**2) / uavg(1)
            tmp1 = asqr + bsqr        
            tmp2 = DSQRT( (asqr - bsqr)**2 + 4.d0 * asqr * btsqr )
            cf_avg   = DSQRT( 0.5d0 * (tmp1 + tmp2) )
            
            vx_avg = uavg(2)/ uavg(1)

            lambdaf_minus_avg = vx_avg - cf_avg 
            lambdaf_plus_avg  = vx_avg + cf_avg


            sL = MIN(lambdaf_minus_L, lambdaf_minus_avg, 0.d0)
            sR = MAX(lambdaf_plus_R, lambdaf_plus_avg, 0.D0)
       
            ! check for max wave speed
            !smax = MAX(smax, ABS(sL), ABS(sR))
            smax = MAX(smax, ABS(vxL)+cfL, ABS(vxR)+cfR, ABS(vx_avg)+cf_avg)
            

            !################### SOME OTHER OPTION FOR CHOOSING THE OUTER WAVE SPEEDS ###########################               

            !sL = MIN(lambdaf_minus_L, lambdaf_minus_avg, 0.d0)
            !sR = MAX(lambdaf_plus_R, lambdaf_plus_avg, 0.D0)

            !sR = MAX(ABS(vx_avg) + cf_avg, ABS(vxR) + cfR , ABS(vxL) + cfL) ! enlarged outer wave-speeds for extra-diffusiveness (see Einfeldt paper) 
            !sL = -sR

            !sL = MIN(lambdaf_minus_L, lambdaf_minus_R)
            !sR = MAX(lambdaf_plus_L, lambdaf_plus_R) 
            
            !sL = MIN(lambdaf_minus_L, lambdaf_minus_R, lambdaf_minus_avg, 0.d0)
            !sR = MAX(lambdaf_plus_L, lambdaf_plus_R, lambdaf_plus_avg, 0.d0)           
            
            !sL = MIN(MIN(lambdaf_minus_L, lambdaf_minus_avg), 0.d0) !, 0.d0)
            !sR = MAX(MAX(lambdaf_plus_R, lambdaf_plus_avg), 0.D0)
           
            !sL = MIN( MIN(lambdaf_minus_L, lambdaf_minus_R) , 0.d0)
            !sR = MAX( MAX(lambdaf_plus_L, lambdaf_plus_R) , 0.d0)       
                        
            !sL = MIN(vxL, vxR, vx_avg) - MIN(cfL, cfR, cf_avg)
            !sR = MAX(vxL, vxR, vx_avg) + MAX(cfL, cfR, cf_avg) 
            
            !####################################################################################################
            
            ! first, check speeds at Riemann fan boundary
            IF(sL .GE. 0.d0) THEN 
                Flux(i,j,:) = FL(:)
                CYCLE
            END IF   

            IF(sR .LE. 0.d0) THEN
                Flux(i,j,:) = FR(:)       
                CYCLE
            END IF            
                               
            ! now compute HLL state vector and flux vector            
            isRL = 1.d0 / (sR - sL)
            uHLL(:) = (sR * uR(:) - sL * uL(:) - FR(:) + FL(:) ) * isRL
            FHLL(:)  = (sR * FL(:) - sL * FR(:) + sR * sL* (uR(:) - uL(:)) ) * isRL        

            ! compute star region continuous variables (i.e. density, momentum and their corresponding fluxes)
            rho_s  = uHLL(1)                      
            momx_s = uHLL(2)
            vx_s   = FHLL(1) / rho_s 
                                  
            irho_s  = 1.d0 / rho_s
            isrho_s = DSQRT(irho_s)
                
            ! compute star region alfven wave speeds
            sL_s = vx_s - ABS(Bx) * isrho_s     
            sR_s = vx_s + ABS(Bx) * isrho_s    

            
            !IF(sL_s .LT. sL .OR. sR_s .GT. sR) THEN
            !    PRINT*,''
            !    PRINT*,'Bad HLLD wave speeds'
            !    PRINT*,'sL, sL_s, vx_s, sR_s, sR = ',sL, sL_s, vx_s, sR_s, sR 
            !    PRINT*,''
            !    PRINT*,''
            !    STOP
            !END IF
            
            !####################################
            ! check for and treat degenrate cases (does not work correctly)  
            !#################################### 
            ! version 1 (use this for better accuracy) 
            !IF(ABS(Bx) .LT. ZERO) THEN
            !    Flux(i,j,:) = FHLL(:)  
            !    CYCLE
            !END IF
            
            ! version 2 (more diffusive than version 1)
            !tmp1 = sR - sL
            !IF(((sL-sL_s) .LT. 1.d-5*tmp1) .OR. ((sR-sR_s) .GT. -1.d-5*tmp1)) THEN
            !    Flux(i,j,:) = FHLL(:)  
            !    CYCLE
            !END IF
            !####################################
            !####################################
                   
            ! compute star region left state
            uL_s(1) = rho_s
            uL_s(2) = momx_s
        
            IF((ABS(sL - sL_s) .LT. ZERO) .OR. (ABS(sL - sR_s) .LT. ZERO)) THEN  ! check for degenerate case (i.e. star region alfven wave speed approaches close to outermost signal speeds)   
                uL_s(3:7) = uL(3:7) ! if degenrate, force the transverse variables to be continuous
            ELSE
                isRL = 1.d0 / ((sL-sL_s) * (sL-sR_s)) 
                tmp1 = Bx * (vx_s - vxL) * isRL
                tmp2 = ( uL(1) * (sL - vxL)**2 - Bx**2 ) * irho_s * isRL
                    
                uL_s(3) = rho_s * vyL -  tmp1 * uL(6)
                uL_s(4) = rho_s * vzL -  tmp1 * uL(7)
                uL_s(5) = uL(5)
                uL_s(6) = tmp2 * uL(6)
                uL_s(7) = tmp2 * uL(7)            
            END IF
           
            ! compute star region right state
            uR_s(1) = rho_s
            uR_s(2) = momx_s
        
            IF((ABS(sR - sL_s) .LT. ZERO) .OR. (ABS(sR - sR_s) .LT. ZERO)) THEN  ! check for degenerate case (i.e. star region alfven wave speed approaches close to outermost signal speeds)   
                uR_s(3:7) = uR(3:7) ! if degenrate, force the transverse variables to be continuous
            ELSE
                isRL = 1.d0 / ((sR-sL_s) * (sR-sR_s)) 
                tmp1 = Bx * (vx_s - vxR) * isRL
                tmp2 = ( uR(1) * (sR - vxR)**2 - Bx**2 ) * irho_s * isRL
                    
                uR_s(3) = rho_s * vyR -  tmp1 * uR(6)
                uR_s(4) = rho_s * vzR -  tmp1 * uR(7)
                uR_s(5) = uR(5)
                uR_s(6) = tmp2 * uR(6)
                uR_s(7) = tmp2 * uR(7)
            END IF
         
            ! compute star region center state
            X  = DSQRT(rho_s) * SIGN(1.d0, Bx)
            iX = 1.d0 / X
                
            uC_s(1) = rho_s
            uC_s(2) = momx_s                                
            uC_s(3) = 0.5d0 * ( uL_s(3) + uR_s(3) + X*( uR_s(6)-uL_s(6) ) )
            uC_s(4) = 0.5d0 * ( uL_s(4) + uR_s(4) + X*( uR_s(7)-uL_s(7) ) )
            uC_s(5) = uL(5)
            uC_s(6) = 0.5d0 * ( uL_s(6) + uR_s(6) + iX*( uR_s(3)-uL_s(3) ) )
            uC_s(7) = 0.5d0 * ( uL_s(7) + uR_s(7) + iX*( uR_s(4)-uL_s(4) ) )
               

            ! finally, compute the HLLD interface flux
            IF(sL_s .GE. 0.d0) THEN
                Flux(i,j,:) = FL(:) + SL * (uL_s(:) - uL(:))
            ELSE IF(sR_s .LE. 0.d0) THEN
                Flux(i,j,:) = FR(:) + SR * (uR_s(:) - uR(:))
            ELSE
                Flux(i,j,1) = FHLL(1)
                Flux(i,j,2) = FHLL(2)
                Flux(i,j,3) = uC_s(3) * vx_s - Bx * uC_s(6)
                Flux(i,j,4) = uC_s(4) * vx_s - Bx * uC_s(7)
                Flux(i,j,5) = 0.d0
                Flux(i,j,6) = vx_s * uC_s(6) - (uC_s(3) * Bx) / rho_s
                Flux(i,j,7) = vx_s * uC_s(7) - (uC_s(4) * Bx) / rho_s   
            END IF           

        END DO    
    END DO   
   
    
    
END SUBROUTINE compute_HLLD_fluxes_1d


! this subroutine checks state vector array for negative densities 
SUBROUTINE final_check_unphysical(ilow, ihi, jlow, jhi, klow, khi)

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi, klow, khi
    INTEGER :: i, j, k

    DO k = klow, khi
        DO j = jlow, jhi
            DO i = ilow, ihi
    
                ! replace negative densities with user-defined floor value
                IF(q_3d(i,j,k,1) .LT. MIN_DENS) THEN
                    PRINT*,'Negative density detected in cell# ',i+my_coord(1)*nx, j+my_coord(2)*ny, k+my_coord(3)*nz
                    q_3d(i,j,k,1) = MIN_DENS 
                END IF
                    
            END DO
        END DO
    END DO

END SUBROUTINE final_check_unphysical


! compute maximum speed (need this for calculating the time step size)
SUBROUTINE get_max_speed_vanleer(max_speed)

    REAL*8, INTENT(INOUT) :: max_speed    
    INTEGER :: i, j, k , imax, jmax, kmax
    REAL*8 :: rho, isrho, vx, vy, vz, bx, by, bz, bsqr, asqr, btsqr, cfx, cfy, cfz, tmp
    
    max_speed = 0.d0

    asqr  = sound_speed**2
            
    ! find the maximum fast mode speed on the grid (the fast mode speed at a cell-interface
    ! is computed using the roe-averaged state at that interface)
    DO k = 1, nz    
      DO j = 1, ny    
        DO i = 1, nx ! no need to include boundary cells for this
        
            rho = q_3d(i,j,k,1) 
            
            isrho = 1.D0 / DSQRT(rho)
            bx  = q_3d(i,j,k,5) * isrho
            by  = q_3d(i,j,k,6) * isrho
            bz  = q_3d(i,j,k,7) * isrho
            btsqr = by**2 + bz**2 
            bsqr  = bx**2 + btsqr 
            tmp = DSQRT( (asqr - bsqr)**2 + 4.d0 * asqr * btsqr )
            
            cfx = DSQRT( 0.5d0 * (asqr + bsqr + tmp) ) ! fast mode velocity x

            
            btsqr = bx**2 + bz**2 
            bsqr  = by**2 + btsqr 
            tmp = DSQRT( (asqr - bsqr)**2 + 4.d0 * asqr * btsqr )
           
            cfy = DSQRT( 0.5d0 * (asqr + bsqr + tmp) ) ! fast mode velocity y

            btsqr = bx**2 + by**2 
            bsqr  = bz**2 + btsqr 
            tmp = DSQRT( (asqr - bsqr)**2 + 4.d0 * asqr * btsqr )
            
            cfz = DSQRT( 0.5d0 * (asqr + bsqr + tmp) ) ! fast mode velocity z
            
            vx = q_3d(i,j,k,2) / rho 
            vy = q_3d(i,j,k,3) / rho
            vz = q_3d(i,j,k,4) / rho

            max_speed = MAX(max_speed, ABS(vx) + cfx, ABS(vy) + cfy, ABS(vz) + cfz)
            
            !IF(max_speed .GT. tmp) THEN
            !    imax = i  
            !    jmax = j
            !    kmax = k
            !END IF    
            
        END DO    
      END DO
    END DO

    max_speed = MAX(max_speed, smax + delv_max)


    !PRINT*,'myrank, maxspeed, Maxspeed cell# = ',myrank, max_speed, imax+my_coord(1)*nx, jmax+my_coord(2)*ny, kmax+my_coord(3)*nz

    IF(max_speed .LT. ZERO) THEN
        PRINT*,'ERROR!! ZERO MAX SPEED ENCOUNTERED!'
        STOP
    END IF


END SUBROUTINE get_max_speed_vanleer


END MODULE VanLeer_Solver_mod