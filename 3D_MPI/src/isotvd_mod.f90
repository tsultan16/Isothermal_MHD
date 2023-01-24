! 3D Isothermal MHD 2nd order TVD Solver
! References: Kim 'et al., Ap. J. 514, 506-519, 1998  (K98)
!             Stone 'et al., Ap. J. S., 178, 137-177 2018 (S18)
!             Mignone, J. Comp., 225, 1427-1441, 2007 (M07)
!
!-----------------------------------------
! Outline of 3D CTU+CT Algorithm from S18:
!-----------------------------------------
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
! *****************
! * CTU Algorithm * 
! *****************
!
! The 3D update at each time step is performed simulatenously in both 
! x,y,z directions using a (sort of) predictor-corrector approach.
!
! 1) Predictor Stage: We first perform directionally split 1/2 time step updates on the 
!                     MHD state variables. The update in a given 
!                     direction is done using the fluxes in the transverse directions.
!                     The cell-face B-field undergoes a half time step CT update.
!                     (Note: The dimensionally split 1/2 time updates require extra "multidimensional
!                     source terms" to be included. These terms are proportional to the derivative o
!                     f B-longitudinal along the given direction).
!
! 2) Corrector Stage: The MHD state variables and cell-face magnetic fields obtained from 
!                     the half-update in the predictor stage are then used to compute "corrector"  
!                     TVD fluxes and cell-corner (cell-edges in 3D) electric electric fields.
!                     These "corrector" fluxes and electric fields are then used to fully update 
!                     the MHD state variables without using any dimensional splitting, i.e. the 
!                     x, y and z TVD flux gradients are simultaneously applied (and the CT magnetic 
!                     field update is already dimensionally unsplit). 
!
!                     The 1D TVD cell-interface fluxes are computed using the 2nd order (via flux-limiters)
!                      Roe-type Riemann solver (K98).
!
!
!
! ************************
! * Protection Algorithm *
! ************************
! When enabled, we concurrently compute "HLLD" fluxes. Each time the state variables get updated,
! we check for negative densities. If we detect a cell containing a negative density, we reupdate
! the state variables in that cell and in it's immediate neighbors using the HLLD fluxes at the 
! boundaries of the bad cell. Computing the HLLD fluxes adds an extra ~15% cost to overall TVD solve
! per time-step. Note: Coputing HLLD fluxes is slightly more expensive than the simpler "HLLE" fluxes
! because HLLE (which only includes the Fast modes) has onle one intermediate state in the star-region,  
! whereas HLLD has 3 (beacause it include both Alfven and Fast modes). 
!
!-------------------------------------------------------------------------------------------------------------


MODULE Iso_TVD_Solver_mod

USE constants_mod
USE grid_data_mod

IMPLICIT NONE


!*****************************************************************************************************************************!
!                                                     LOCAL VARIABLES                                                         !
!*****************************************************************************************************************************!


! Solver Parameters
INTEGER, PARAMETER :: nwaves = 6
REAL*8, PARAMETER  :: ZERO = 1.D-20           ! density floor
REAL*8, PARAMETER :: limiter_beta = 1.D0      ! Sweby limiter: 1 <= beta <= 2
REAL*8, PARAMETER :: MIN_DENS = 1.D-20        ! density floor
LOGICAL, PARAMETER :: flux_protection = .FALSE.

! Solver Local variables 
REAL*8, ALLOCATABLE :: q(:,:,:),             & ! state-vector 2D work array
                       qwork_3d(:,:,:,:),    & ! state-vector 3D work array
                       q_3d_prot(:,:,:,:),   & ! storage for state vector protection backup
                       bface_work(:,:,:,:),  & ! cell-face magnetic field work array 
                       Fs(:,:,:,:,:),        & ! time-averaged flux vector
                       Fs_buff(:,:,:,:,:),   & ! flux buffer for storing intermediate/predictor fluxes
                       Fs_prot(:,:,:,:,:),   & ! storage for protection (HLLD) fluxes
                       Fs_buff_prot(:,:,:,:,:), & ! storage for protection (HLLD) fluxes
                       Fs_2d(:,:,:),         & ! flux buffer for rotations
                       Fs_2d_prot(:,:,:),    & ! 2d array for protection flux storage
                       qavg(:,:,:),          & ! Roe-averaged state    
                       Rk(:,:,:,:),          & ! Roe-averaged right eigenvectors
                       ck(:,:,:),            & ! characteristic variables (i.e. projection of state vector onto the left eigenvectors)   
                       eigenvalues(:,:,:),   & ! Roe-averaged eigenvalues
                       wave_speeds(:,:,:),   & ! Roe-averaged wave speeds
                       emf_3d(:,:,:,:,:),    & ! cell-center reference and cell-face electric field z-component
                       emf_3d_prot(:,:,:,:,:), & ! storage protection emf
                       emf_corner_prot(:,:,:,:)  ! cell-corner CT electric field

LOGICAL, ALLOCATABLE :: protection_flags(:,:,:)


REAL*8, PARAMETER :: eps(6) = (/ 0.2, 0.0, 0.2, 0.2, 0.0 ,0.2/) ! wave dissipation constants

REAL*8 :: dt, dx

INTEGER :: nmax


!*****************************************************************************************************************************!
!                                                     SUBROUTINES                                                             !
!*****************************************************************************************************************************!

CONTAINS


SUBROUTINE init_solver_tvd()
    
    ! make sure we have at least 5 boundary cells
    IF(nb .LT. 5) THEN
        PRINT*,'ERROR! Need nb >= 5 for TVD solver.'
        STOP
    END IF
    
    CALL allocate_local_variables()


END SUBROUTINE init_solver_tvd


SUBROUTINE destroy_solver_tvd()

   CALL deallocate_local_variables()
   !CLOSE(UNIT=12)
   
END SUBROUTINE destroy_solver_tvd

! allocate local variables/work arrays
SUBROUTINE allocate_local_variables()

    INTEGER :: total_grid_memory

    nmax = MAX(nx*nranks_x,ny*nranks_y,nz*nranks_z)

    ! set cell size  
    dx = 1.d0 / DBLE(nmax)
    
    ALLOCATE(q(1-nb:nmax+nb,1-nb:nmax+nb,7))
    ALLOCATE(qwork_3d(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,7))
    ALLOCATE(bface_work(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb, 3))
    ALLOCATE(Fs(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,7,3))
    ALLOCATE(Fs_buff(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,7,3))
    ALLOCATE(Fs_2d(1-nb:nmax+nb,1-nb:nmax+nb,7))
    ALLOCATE(qavg(1-nb:nmax+nb,1-nb:nmax+nb,7))
    ALLOCATE(Rk(1-nb:nmax+nb,1-nb:nmax+nb,nwaves, nwaves))
    ALLOCATE(ck(1-nb:nmax+nb,1-nb:nmax+nb,nwaves))
    ALLOCATE(eigenvalues(1-nb:nmax+nb,1-nb:nmax+nb,nwaves))
    ALLOCATE(wave_speeds(1-nb:nmax+nb,1-nb:nmax+nb,3))
    ALLOCATE(emf_3d(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,3,3))

    IF(flux_protection) THEN

        ALLOCATE(q_3d_prot(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,7))
        ALLOCATE(Fs_prot(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,7,3))
        ALLOCATE(Fs_buff_prot(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,7,3))
        ALLOCATE(Fs_2d_prot(1-nb:nmax+nb,1-nb:nmax+nb,7))
        ALLOCATE(emf_3d_prot(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,3,3))
        ALLOCATE(emf_corner_prot(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,3))
        ALLOCATE(protection_flags(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb))

    END IF


    total_grid_memory = SIZEOF(q) + SIZEOF(qwork_3d) + SIZEOF(bface_work) + &
                        SIZEOF(Fs) + SIZEOF(Fs_buff) + SIZEOF(Fs_2d) + &
                        SIZEOF(qavg) + SIZEOF(Rk) + SIZEOF(ck) + SIZEOF(eigenvalues) + &
                        SIZEOF(wave_speeds) + SIZEOF(emf_3d) 
    
    IF(flux_protection) THEN
    
        total_grid_memory =  total_grid_memory + SIZEOF(q_3d_prot) + SIZEOF(Fs_prot) + &
                             SIZEOF(Fs_buff_prot) + SIZEOF(Fs_2d_prot) + &
                             SIZEOF(emf_3d_prot) + SIZEOF(emf_corner_prot)
    END IF

    IF(myrank .EQ. 0) THEN
    
    PRINT*,''
    PRINT*,'#####################################################################'
    PRINT*,'Total Memory for work array storage (Mb) = ',total_grid_memory*1e-6                        
    PRINT*,'#####################################################################'
    PRINT*,''

    END IF
    

END SUBROUTINE allocate_local_variables


SUBROUTINE deallocate_local_variables()


    DEALLOCATE(q, qwork_3d, bface_work)
    DEALLOCATE(Fs, Fs_buff, qavg, Rk, ck, eigenvalues, wave_speeds)
    DEALLOCATE(emf_3d)
    DEALLOCATE(Fs_2d)
    
    IF(flux_protection) THEN 
        DEALLOCATE(q_3d_prot, Fs_buff_prot, Fs_prot, Fs_2d_prot)    
        DEALLOCATE(emf_3d_prot, emf_corner_prot)
        DEALLOCATE(protection_flags)
    END IF
    
    
END SUBROUTINE deallocate_local_variables


!*****************************************************************************************************************************!

! Top-level Isothermal TVD MHD solve routine (using CTU + CT scheme)
SUBROUTINE tvd_solve_3d(dt_in)

    REAL*8, INTENT(INOUT) :: dt_in
    INTEGER :: i, j, offset  
    
    
    dt = dt_in
    
    ! clear flux arrays (maybe unnecessary)
    Fs = 0.d0
    Fs_buff = 0.d0
    Fs_prot = 0.d0

    ! clear electric field array (may not be necessary)
    emf_3d = 0.d0
    emf_corner = 0.d0  

    ! First save initial state in work arrays
    qwork_3d = q_3d
    bface_work = bface_3d

    offset = 0
    
    !***********************************************************************************
    ! Step 1: Compute 1/2 dt TVD fluxes along x, y and z directions using 1D TVD solver
    !***********************************************************************************
    
    IF(print_debug) PRINT*,'Computing predictor TVD x-fluxes.'
    
    CALL compute_xfluxes_3d(1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset, .FALSE.)  ! send false flag since this is the predictor stage

    IF(print_debug) PRINT*,'Computing predictor TVD y-fluxes.'
    
    CALL compute_yfluxes_3d(1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset, .FALSE.)  ! send false flag since this is the predictor stage

    IF(print_debug) PRINT*,'Computing predictor TVD z-fluxes.'
    
    CALL compute_zfluxes_3d(1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset, .FALSE.)  ! send false flag since this is the predictor stage
    

    !*******************************************************************************
    ! Step 2: Compute cell-centered reference and cell-corner electric fields  
    !*******************************************************************************    
    
    IF(print_debug) PRINT*,'Computing predictor CT cell center reference electric field.'

    CALL compute_cell_center_emf_3d(.FALSE.)  ! send false flag since this is the predictor stage

    IF(print_debug) PRINT*,'Computing predictor CT cell corner electric field.'

    CALL compute_cell_corner_emf_3d(emf_3d, emf_corner, 1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset)
    IF(flux_protection) CALL compute_cell_corner_emf_3d(emf_3d_prot, emf_corner_prot, 1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset)

    !*************************************************************************************
    ! Step 3: Apply 1/2 dt CT update to cell-face magnetic fields and compute corrector
    !         cell-center reference electric field
    !************************************************************************************
    
    IF(print_debug) PRINT*,'Computing CT half update.'

    CALL update_CT_3d(1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset, .FALSE.) ! send false flag since this is the predictor stage
    
    IF(print_debug) PRINT*,'Computing corrector CT cell-center reference electric field.' 

    CALL compute_cell_center_emf_3d(.TRUE.)   ! send true flag since this is the corrector stage
    
    
    !*******************************************************************************
    ! Step 4: Apply dimensionally split 1/2 dt updates to state variables using
    !         transverse fluxes (including multidimensional source terms) immediately 
    !         followed by the "corrector" flux calculation.
    !*******************************************************************************
    
    offset = 2
    
    IF(print_debug) PRINT*,'Computing state variables half update and final TVD fluxes in x direction.'

    ! 1/2 dt update along x
    CALL halfupdate_x_3d(Fs_buff)
     
    ! compute corrector TVD x-fluxes     
    CALL compute_xfluxes_3d(1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset, .TRUE.)     ! send true flag since this is the corrector stage
    

    IF(print_debug) PRINT*,'Computing state variables half update and final TVD fluxes in y direction.'
    
    ! Restore intial state vector
    q_3d = qwork_3d

    ! 1/2 dt update along y
    CALL halfupdate_y_3d(Fs_buff)
     
    ! compute corrector TVD y-fluxes     
    CALL compute_yfluxes_3d(1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset, .TRUE.)     ! send true flag since this is the corrector stage    


    IF(print_debug) PRINT*,'Computing state variables half update and final TVD fluxes in z direction.'

    ! Restore intial state vector
    q_3d = qwork_3d
    
    ! 1/2 dt update along z    
    CALL halfupdate_z_3d(Fs_buff)
  
    ! compute corrector TVD z-fluxes
    CALL compute_zfluxes_3d(1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset, .TRUE.)     ! send true flag since this is the corrector stage


    !*********************************************************************************
    ! Step 5: Compute the final cell-corner electric fields   
    !*********************************************************************************
    
    ! Restore intial state vector
    q_3d = qwork_3d
    
    IF(print_debug) PRINT*,'Computing corrector CT cell-corner electric field.' 

    CALL compute_cell_corner_emf_3d(emf_3d, emf_corner, 1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset)
    IF(flux_protection) CALL compute_cell_corner_emf_3d(emf_3d_prot, emf_corner_prot, 1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset)
    
    !*******************************************************************************
    ! Step 6: Apply full dt update to density, momentum and cell face magnetic 
    !         fields using the new fluxes and cell-corner electic fields
    !******************************************************************************* 
    
    IF(print_debug) PRINT*,'Computing final unsplit update of state variables and CT.'
        
    CALL update_state_vector_3d(Fs)
   
    ! restore original cell face magnetic fields 
    bface_3d = bface_work
        
    CALL update_CT_3d(1+offset, nx-offset, 1+offset, ny-offset, 1+offset, nz-offset, .TRUE.) ! send true flag since this is the corrector stage
    
    IF(print_debug) PRINT*,'Update for time step completed.'
    
 
    
END SUBROUTINE tvd_solve_3d


! Computes 1/2 dt (precictor) TVD x-fluxes
SUBROUTINE compute_xfluxes_3d(ilow, ihi, jlow, jhi, klow, khi, corrector_stage) 

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi, klow, khi
    LOGICAL, INTENT(IN) :: corrector_stage
    INTEGER :: i,j,k
    REAL*8 :: dtdx


    dtdx = dt / dx   ! for half-step

    IF(.NOT. corrector_stage)  dtdx = 0.5d0 * dtdx ! for half-update during predictor stage
 
 
    ! Loop over z planes    
    DO k = klow-nb, khi+nb
 
        IF(print_debug) PRINT*,' Z plane index,  k = ', k
        
        ! first, copy state vector into our work array 
        DO j = jlow-nb, jhi+nb 
            DO i = ilow-nb, ihi+nb        
                q(i,j,1) = q_3d(i,j,k,1)
                q(i,j,2) = q_3d(i,j,k,2)
                q(i,j,3) = q_3d(i,j,k,3)
                q(i,j,4) = q_3d(i,j,k,4)
                q(i,j,5) = bface_3d(i,j,k,1) !q_3d(i,j,k,5)
                q(i,j,6) = q_3d(i,j,k,6)
                q(i,j,7) = q_3d(i,j,k,7)
            END DO
        END DO

       
        IF(riemann_solver_type .EQ. 1) THEN
        
            IF(print_debug) PRINT*,'Computing HLLD fluxes.' 

            CALL compute_HLLE_fluxes_1d(Fs_2d, ilow, ihi, jlow, jhi)
        
        ELSE IF(riemann_solver_type .EQ. 2) THEN
 
 
            IF(print_debug) PRINT*,'Computing Roe average state.' 
    
            ! now sweep through strips along x-direction and compute x-interface fluxes   
            CALL compute_roe_avg_state_1d(ilow, ihi, jlow, jhi)
    
            IF(print_debug) PRINT*,'Computing eigenvalues.' 
    
            CALL compute_eigenvalues_1d(ilow, ihi, jlow, jhi)
            
            IF(print_debug) PRINT*,'Computing eigenvectors.'
        
            CALL compute_eigenvectors_1d(ilow, ihi, jlow, jhi)

            IF(print_debug) PRINT*,'Computing x-fluxes.'

            CALL compute_TVD_fluxes_1d(Fs_2d, ilow, ihi, jlow, jhi, dtdx)

            IF(flux_protection) CALL compute_HLLE_fluxes_1d(Fs_2d_prot, ilow, ihi, jlow, jhi)

        END IF
                
        ! for predictor stage, store TVD fluxes in the buffer array
        IF(.NOT. corrector_stage) THEN

            DO j = jlow-nb, jhi+nb 
                DO i = ilow-nb, ihi+nb 
                    Fs_buff(i,j,k,:,1) = Fs_2d(i,j,:)
                END DO
            END DO
        
            ! this is a good place to store the cell-face transverse electric fields
            DO j = jlow-nb, jhi+nb
                DO i = ilow-nb+1, ihi+nb-2     
                    emf_3d(i,j,k,3,1) = -Fs_buff(i,j,k,6,1)   ! Flux-x(By)_i+1/2,j,k = -Ez_i+1/2,j,k 
                    emf_3d(i,j,k,2,1) = Fs_buff(i,j,k,7,1)    ! Flux-x(Bz)_i+1/2,j,k = Ey_i+1/2,j,k
                    
                END DO
            END DO    
    
            IF(flux_protection) THEN
    
                DO j = jlow-nb, jhi+nb 
                    DO i = ilow-nb, ihi+nb 
                        Fs_buff_prot(i,j,k,:,1) = Fs_2d_prot(i,j,:)
                    END DO
                END DO
        
                ! this is a good place to store the cell-face transverse electric fields
                DO j = jlow-nb, jhi+nb
                    DO i = ilow-nb+1, ihi+nb-2     
                        emf_3d_prot(i,j,k,3,1) = -Fs_buff_prot(i,j,k,6,1)   ! Flux-x(By)_i+1/2,j,k = -Ez_i+1/2,j,k 
                        emf_3d_prot(i,j,k,2,1) = Fs_buff_prot(i,j,k,7,1)    ! Flux-x(Bz)_i+1/2,j,k = Ey_i+1/2,j,k
                    
                    END DO
                END DO    
    
            END IF
    
        ELSE
        
            DO j = jlow-nb, jhi+nb 
                DO i = ilow-nb, ihi+nb 
                    Fs(i,j,k,:,1) = Fs_2d(i,j,:)
                END DO
            END DO
            
            ! this is a good place to store the cell-face transverse electric fields
            DO j = jlow-nb, jhi+nb
                DO i = ilow-nb+1, ihi+nb-2   
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
                    DO i = ilow-nb+1, ihi+nb-2     
                        emf_3d_prot(i,j,k,3,1) = -Fs_prot(i,j,k,6,1)   ! Flux-x(By)_i+1/2,j,k = -Ez_i+1/2,j,k 
                        emf_3d_prot(i,j,k,2,1) = Fs_prot(i,j,k,7,1)    ! Flux-x(Bz)_i+1/2,j,k = Ey_i+1/2,j,k
                    
                    END DO
                END DO    
    
            END IF
    
        END IF
    
    END DO
    

END SUBROUTINE compute_xfluxes_3d 



! Computes 1/2 dt (predictor) TVD y-fluxes
SUBROUTINE compute_yfluxes_3d(ilow, ihi, jlow, jhi, klow, khi, corrector_stage) 

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi, klow, khi
    LOGICAL, INTENT(IN) :: corrector_stage
    INTEGER :: i,j,k
    REAL*8 :: dtdx


    dtdx = dt / dx      
    
    IF(.NOT. corrector_stage)  dtdx = 0.5d0 * dtdx ! for half-update during predictor stage

    ! loop over z planes
    DO k = klow-nb, khi+nb
    
        IF(print_debug) PRINT*,' Z plane index,  k = ', k
            
        ! first rotate the stave vector array so that fastest index runs along y-direction
        ! (also need to permute  x->y->z components for momentum and B field)
        DO i = ilow-nb, ihi+nb 
            DO j = jlow-nb, jhi+nb
        
                q(j,i,1) = q_3d(i,j,k,1)
                q(j,i,2) = q_3d(i,j,k,3)
                q(j,i,3) = q_3d(i,j,k,4)
                q(j,i,4) = q_3d(i,j,k,2)
                q(j,i,5) = bface_3d(i,j,k,2) ! q_3d(i,j,k,6)
                q(j,i,6) = q_3d(i,j,k,7)
                q(j,i,7) = q_3d(i,j,k,5)
        
            END DO
        END DO
     
        
        IF(riemann_solver_type .EQ. 1) THEN
        
            IF(print_debug) PRINT*,'Computing HLLD fluxes.' 

            CALL compute_HLLE_fluxes_1d(Fs_2d, jlow, jhi, ilow, ihi)
        
        ELSE IF(riemann_solver_type .EQ. 2) THEN
        
            IF(print_debug) PRINT*,'Computing Roe average state.' 

            ! now sweep through strips along y-direction and compute y-interface fluxes   
            CALL compute_roe_avg_state_1d(jlow, jhi, ilow, ihi)
    
            IF(print_debug) PRINT*,'Computing eigenvalues.'
     
            CALL compute_eigenvalues_1d(jlow, jhi, ilow, ihi)
       
            IF(print_debug) PRINT*,'Computing eigenvectors.'
    
            CALL compute_eigenvectors_1d(jlow, jhi, ilow, ihi)

            IF(print_debug) PRINT*,'Computing y-fluxes.'

            CALL compute_TVD_fluxes_1d(Fs_2d,jlow, jhi, ilow, ihi, dtdx)

            IF(flux_protection) CALL compute_HLLE_fluxes_1d(Fs_2d_prot, jlow, jhi, ilow, ihi)

        END IF

        ! for predictor stage, store TVD fluxes in the buffer array
        IF(.NOT. corrector_stage) THEN
            
            ! Rotate the y-fluxes in flux buffer and store them in main flux array 
            DO j = jlow-nb, jhi+nb 
                DO i = ilow-nb, ihi+nb
        
                    Fs_buff(i,j,k,1,2) = Fs_2d(j,i,1)
                    Fs_buff(i,j,k,3,2) = Fs_2d(j,i,2)
                    Fs_buff(i,j,k,4,2) = Fs_2d(j,i,3)
                    Fs_buff(i,j,k,2,2) = Fs_2d(j,i,4)
                    Fs_buff(i,j,k,6,2) = Fs_2d(j,i,5)
                    Fs_buff(i,j,k,7,2) = Fs_2d(j,i,6)         
                    Fs_buff(i,j,k,5,2) = Fs_2d(j,i,7)         
            
                END DO
            END DO
    
            ! this is a good place to store the cell-face transverse electric fields
            DO j = jlow-nb+1, jhi+nb-2
                DO i = ilow-nb, ihi+nb   
                    emf_3d(i,j,k,1,1) = -Fs_buff(i,j,k,7,2)   ! Flux-y(Bz)_i,j+1/2,k = -Ex_i,j+1/2,k
                    emf_3d(i,j,k,3,2) = Fs_buff(i,j,k,5,2)  ! Flux-y(Bx)_i,j+1/2,k = -Ez_i,j+1/2,k
                END DO
            END DO      
            
            IF(flux_protection) THEN
 
                ! Rotate the y-fluxes in flux buffer and store them in main flux array 
                DO j = jlow-nb, jhi+nb 
                    DO i = ilow-nb, ihi+nb
        
                        Fs_buff_prot(i,j,k,1,2) = Fs_2d_prot(j,i,1)
                        Fs_buff_prot(i,j,k,3,2) = Fs_2d_prot(j,i,2)
                        Fs_buff_prot(i,j,k,4,2) = Fs_2d_prot(j,i,3)
                        Fs_buff_prot(i,j,k,2,2) = Fs_2d_prot(j,i,4)
                        Fs_buff_prot(i,j,k,6,2) = Fs_2d_prot(j,i,5)
                        Fs_buff_prot(i,j,k,7,2) = Fs_2d_prot(j,i,6)         
                        Fs_buff_prot(i,j,k,5,2) = Fs_2d_prot(j,i,7)         
            
                    END DO
                END DO
    
                ! this is a good place to store the cell-face transverse electric fields
                DO j = jlow-nb+1, jhi+nb-2
                    DO i = ilow-nb, ihi+nb   
                        emf_3d_prot(i,j,k,1,1) = -Fs_buff_prot(i,j,k,7,2)   ! Flux-y(Bz)_i,j+1/2,k = -Ex_i,j+1/2,k
                        emf_3d_prot(i,j,k,3,2) = Fs_buff_prot(i,j,k,5,2)  ! Flux-y(Bx)_i,j+1/2,k = -Ez_i,j+1/2,k
                    END DO
                END DO      
 
            END IF

        ELSE 

            ! Rotate the y-fluxes in flux buffer and store them in main flux array 
            DO j = jlow-nb, jhi+nb 
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
            DO j = jlow-nb+1, jhi+nb-2
                DO i = ilow-nb, ihi+nb  
                    emf_3d(i,j,k,1,1) = -Fs(i,j,k,7,2)   ! Flux-y(Bz)_i,j+1/2,k = -Ex_i,j+1/2,k
                    emf_3d(i,j,k,3,2) = Fs(i,j,k,5,2)  ! Flux-y(Bx)_i,j+1/2,k = -Ez_i,j+1/2,k
                END DO
            END DO      
            
            IF(flux_protection) THEN
 
                ! Rotate the y-fluxes in flux buffer and store them in main flux array 
                DO j = jlow-nb, jhi+nb 
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
                DO j = jlow-nb+1, jhi+nb-2
                    DO i = ilow-nb, ihi+nb  
                        emf_3d_prot(i,j,k,1,1) = -Fs_prot(i,j,k,7,2)   ! Flux-y(Bz)_i,j+1/2,k = -Ex_i,j+1/2,k
                        emf_3d_prot(i,j,k,3,2) = Fs_prot(i,j,k,5,2)  ! Flux-y(Bx)_i,j+1/2,k = -Ez_i,j+1/2,k
                    END DO
                END DO                      
 
            END IF
    
        END IF

    END DO

END SUBROUTINE compute_yfluxes_3d 



! Computes 1/2 dt (predictor) TVD z-fluxes
SUBROUTINE compute_zfluxes_3d(ilow, ihi, jlow, jhi, klow, khi, corrector_stage) 

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi, klow, khi
    LOGICAL, INTENT(IN) :: corrector_stage
    INTEGER :: i,j,k
    REAL*8 :: dtdx


    dtdx = dt / dx      

    IF(.NOT. corrector_stage)  dtdx = 0.5d0 * dtdx ! for half-update during predictor stage

    
    ! loop over y planes
    DO j = jlow-nb, jhi+nb  
    
        IF(print_debug) PRINT*,' Y plane index:  j = ', j
    
        ! first rotate the stave vector array so that fastest index runs along z-direction 
        ! (also need to permute x->y->z components for momentum and B field)
        DO i = ilow-nb, ihi+nb 
            DO k = klow-nb, khi+nb
        
                q(k,i,1) = q_3d(i,j,k,1)
                q(k,i,2) = q_3d(i,j,k,4)
                q(k,i,3) = q_3d(i,j,k,2)
                q(k,i,4) = q_3d(i,j,k,3)
                q(k,i,5) = bface_3d(i,j,k,3)! q_3d(i,j,k,7)
                q(k,i,6) = q_3d(i,j,k,5)
                q(k,i,7) = q_3d(i,j,k,6)
         
            END DO
        END DO
        
        
        
        IF(riemann_solver_type .EQ. 1) THEN
        
            CALL compute_HLLE_fluxes_1d(Fs_2d,  klow, khi, ilow, ihi)
        
        ELSE IF(riemann_solver_type .EQ. 2) THEN
        
            IF(print_debug) PRINT*,'Computing Roe average state.' 

            ! now sweep through strips along y-direction and compute y-interface fluxes   
            CALL compute_roe_avg_state_1d(klow, khi, ilow, ihi)
    
            IF(print_debug) PRINT*,'Computing eigenvalues.'
     
            CALL compute_eigenvalues_1d(klow, khi, ilow, ihi)
       
            IF(print_debug) PRINT*,'Computing eigenvectors.'
    
            CALL compute_eigenvectors_1d(klow, khi, ilow, ihi)

            IF(print_debug) PRINT*,'Computing z-fluxes.'

            CALL compute_TVD_fluxes_1d(Fs_2d, klow, khi, ilow, ihi, dtdx)

            IF(flux_protection) CALL compute_HLLE_fluxes_1d(Fs_2d_prot,  klow, khi, ilow, ihi)

        END IF

        ! for predictor stage, store TVD fluxes in the buffer array
        IF(.NOT. corrector_stage) THEN
        
            ! Rotate the z-fluxes in flux buffer and store them in main flux array 
            DO k = klow-nb, khi+nb 
                DO i = ilow-nb, ihi+nb
        
                    Fs_buff(i,j,k,1,3) = Fs_2d(k,i,1)
                    Fs_buff(i,j,k,4,3) = Fs_2d(k,i,2)
                    Fs_buff(i,j,k,2,3) = Fs_2d(k,i,3)
                    Fs_buff(i,j,k,3,3) = Fs_2d(k,i,4)
                    Fs_buff(i,j,k,7,3) = Fs_2d(k,i,5)
                    Fs_buff(i,j,k,5,3) = Fs_2d(k,i,6)         
                    Fs_buff(i,j,k,6,3) = Fs_2d(k,i,7)         
            
                END DO
            END DO
    
            ! this is a good place to store the cell-face transverse electric fields
            DO k = klow-nb+1, khi+nb-2
                DO i = ilow-nb, ihi+nb   
                    emf_3d(i,j,k,2,2) = -Fs_buff(i,j,k,5,3)   ! Flux-z(Bx)_i,j,k+1/2 = Ey_i,j,k+1/2
                    emf_3d(i,j,k,1,2) = Fs_buff(i,j,k,6,3)    ! Flux-z(By)_i,j,k+1/2 = Ex_i,j,k+1/2              
                END DO
            END DO      

            IF(flux_protection) THEN
            
                ! Rotate the z-fluxes in flux buffer and store them in main flux array 
                DO k = klow-nb, khi+nb 
                    DO i = ilow-nb, ihi+nb
        
                        Fs_buff_prot(i,j,k,1,3) = Fs_2d_prot(k,i,1)
                        Fs_buff_prot(i,j,k,4,3) = Fs_2d_prot(k,i,2)
                        Fs_buff_prot(i,j,k,2,3) = Fs_2d_prot(k,i,3)
                        Fs_buff_prot(i,j,k,3,3) = Fs_2d_prot(k,i,4)
                        Fs_buff_prot(i,j,k,7,3) = Fs_2d_prot(k,i,5)
                        Fs_buff_prot(i,j,k,5,3) = Fs_2d_prot(k,i,6)         
                        Fs_buff_prot(i,j,k,6,3) = Fs_2d_prot(k,i,7)         
            
                    END DO
                END DO
    
                ! this is a good place to store the cell-face transverse electric fields
                DO k = klow-nb+1, khi+nb-2
                    DO i = ilow-nb, ihi+nb   
                        emf_3d_prot(i,j,k,2,2) = -Fs_buff_prot(i,j,k,5,3)   ! Flux-z(Bx)_i,j,k+1/2 = Ey_i,j,k+1/2
                        emf_3d_prot(i,j,k,1,2) = Fs_buff_prot(i,j,k,6,3)    ! Flux-z(By)_i,j,k+1/2 = Ex_i,j,k+1/2              
                    END DO
                END DO      

            END IF

        ELSE 
    
        
            ! Rotate the z-fluxes in flux buffer and store them in main flux array 
            DO k = klow-nb, khi+nb 
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
            DO k = klow-nb+1, khi+nb-2
                DO i = ilow-nb, ihi+nb   
                    emf_3d(i,j,k,2,2) = -Fs(i,j,k,5,3)   ! Flux-z(Bx)_i,j,k+1/2 = Ey_i,j,k+1/2
                    emf_3d(i,j,k,1,2) = Fs(i,j,k,6,3)    ! Flux-z(By)_i,j,k+1/2 = Ex_i,j,k+1/2              
                END DO
            END DO      
            
            IF(flux_protection) THEN

                ! Rotate the z-fluxes in flux buffer and store them in main flux array 
                DO k = klow-nb, khi+nb 
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
                DO k = klow-nb+1, khi+nb-2
                    DO i = ilow-nb, ihi+nb   
                        emf_3d_prot(i,j,k,2,2) = -Fs_prot(i,j,k,5,3)   ! Flux-z(Bx)_i,j,k+1/2 = Ey_i,j,k+1/2
                        emf_3d_prot(i,j,k,1,2) = Fs_prot(i,j,k,6,3)    ! Flux-z(By)_i,j,k+1/2 = Ex_i,j,k+1/2              
                    END DO
                END DO    

            END IF
    
        END IF

    END DO

END SUBROUTINE compute_zfluxes_3d 


! Compute the cell-centered refernece emf: Ez_i,j
SUBROUTINE compute_cell_center_emf_3d(corrector_stage)

    LOGICAL, INTENT(IN) :: corrector_stage
    INTEGER :: i, j, k
    REAL*8 :: rho, vx, vy, vz, momx, momy, momz, Bx, By, Bz, dtdx
    
    
    ! For the corrector stage, wee need the half-updated cell-center momentum
    IF(corrector_stage) THEN
            
            
        dtdx = 0.5d0 * dt / dx
        
        DO k = 1-nb+3, nz+nb-3
            DO j = 1-nb+3, ny+nb-3
                DO i = 1-nb+3, nx+nb-3     

                    ! To get cell-centered emf, we need cell-centered velocity at half-time step.
                    ! Can get those by applying an unsplit conservative half-update of density and 
                    ! momentum using the already available 1/2 dt fluxes. 
                    rho = q_3d(i,j,k,1) - dtdx * (Fs_buff(i,j,k,1,1) - Fs_buff(i-1,j,k,1,1) +  &
                          Fs_buff(i,j,k,1,2) - Fs_buff(i,j-1,k,1,2) + Fs_buff(i,j,k,1,3) - Fs_buff(i,j,k-1,1,3))  
                    
                    momx = q_3d(i,j,k,2) - dtdx * (Fs_buff(i,j,k,2,1) - Fs_buff(i-1,j,k,2,1) + &
                           Fs_buff(i,j,k,2,2) - Fs_buff(i,j-1,k,2,2) + Fs_buff(i,j,k,2,3) - Fs_buff(i,j,k-1,2,3) )  
               
                    momy = q_3d(i,j,k,3) - dtdx * (Fs_buff(i,j,k,3,1) - Fs_buff(i-1,j,k,3,1) + &
                           Fs_buff(i,j,k,3,2) - Fs_buff(i,j-1,k,3,2) + Fs_buff(i,j,k,3,3) - Fs_buff(i,j,k-1,3,3) )

                    momz = q_3d(i,j,k,4) - dtdx * (Fs_buff(i,j,k,4,1) - Fs_buff(i-1,j,k,4,1) + &
                           Fs_buff(i,j,k,4,2) - Fs_buff(i,j-1,k,4,2) + Fs_buff(i,j,k,4,3) - Fs_buff(i,j,k-1,4,3) )
                       
                    vx = momx / rho
                    vy = momy / rho
                    vz = momz / rho
         
                    ! Get magnetic field components at half-time step by averaging over cell-face
                    ! values obtained from step 3.  
                    Bx = 0.5d0 * (bface_3d(i,j,k,1) + bface_3d(i-1,j,k,1))         
                    By = 0.5d0 * (bface_3d(i,j,k,2) + bface_3d(i,j-1,k,2))    
                    Bz = 0.5d0 * (bface_3d(i,j,k,3) + bface_3d(i,j,k-1,3))    

                    emf_3d(i,j,k,1,3) = By * vz - Bz * vy    ! Ex_i,j,k          
                    emf_3d(i,j,k,2,3) = Bz * vx - Bx * vz    ! Ey_i,j,k      
                    emf_3d(i,j,k,3,3) = Bx * vy - By * vx    ! Ez_i,j,k      
        
                END DO
            END DO    
        END DO
        
        
        IF(flux_protection) THEN
        
             DO k = 1-nb+3, nz+nb-3
                DO j = 1-nb+3, ny+nb-3
                    DO i = 1-nb+3, nx+nb-3     

                        emf_3d_prot(i,j,k,1,3) = emf_3d(i,j,k,1,3)
                        emf_3d_prot(i,j,k,2,3) = emf_3d(i,j,k,2,3)
                        emf_3d_prot(i,j,k,3,3) = emf_3d(i,j,k,3,3)
                        

                    END DO
                END DO    
            END DO
        
        END IF
        
        
    ELSE
        
        DO k = 1-nb, nz+nb
            DO j = 1-nb, ny+nb
                DO i = 1-nb, nx+nb        
         
                    emf_3d(i,j,k,1,3) = ( q_3d(i,j,k,6) * q_3d(i,j,k,4)  - q_3d(i,j,k,7) * q_3d(i,j,k,3) ) / q_3d(i,j,k,1)   ! Ex_i,j,k          
                    emf_3d(i,j,k,2,3) = ( q_3d(i,j,k,7) * q_3d(i,j,k,2)  - q_3d(i,j,k,5) * q_3d(i,j,k,4) ) / q_3d(i,j,k,1)   ! Ey_i,j,k      
                    emf_3d(i,j,k,3,3) = ( q_3d(i,j,k,5) * q_3d(i,j,k,3)  - q_3d(i,j,k,6) * q_3d(i,j,k,2) ) / q_3d(i,j,k,1)   ! Ez_i,j,k 
                
                END DO
            END DO    
        END DO


        IF(flux_protection) THEN
        
             DO k = 1-nb, nz+nb
                DO j = 1-nb, ny+nb
                    DO i = 1-nb, nx+nb     

                        emf_3d_prot(i,j,k,1,3) = emf_3d(i,j,k,1,3)
                        emf_3d_prot(i,j,k,2,3) = emf_3d(i,j,k,2,3)
                        emf_3d_prot(i,j,k,3,3) = emf_3d(i,j,k,3,3)
                        

                    END DO
                END DO    
            END DO
        
        END IF
        
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
    DO k = klow-nb+1, khi+nb-2
        DO j = jlow-nb+1, jhi+nb-2
            DO i = ilow-nb+1, ihi+nb-2  
        
                ! get normal velocity at x-interface (roe-average value)
                vx = 0.5d0 * (q_3d(i,j,k,2)/q_3d(i,j,k,1) + q_3d(i+1,j,k,2)/q_3d(i+1,j,k,1) )
            
                ! get normal velocity at y-interface (roe-average value)
                vy = 0.5d0 * (q_3d(i,j,k,3)/q_3d(i,j,k,1) + q_3d(i,j+1,k,3)/q_3d(i,j+1,k,1) )
    
                ! compute upwinded interface gradients of electric field at x-interface
                dyL_14 = idx2 * (emf_3d_in(i,j,k,3,2) - emf_3d_in(i,j,k,3,3)) 
                dyR_14 = idx2 * (emf_3d_in(i+1,j,k,3,2) - emf_3d_in(i+1,j,k,3,3)) 
                dyL_34 = idx2 * (emf_3d_in(i,j+1,k,3,3) - emf_3d_in(i,j,k,3,2)) 
                dyR_34 = idx2 * (emf_3d_in(i+1,j+1,k,3,3) - emf_3d_in(i+1,j,k,3,2))

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
                dxL_14 = idx2 * (emf_3d_in(i,j,k,3,1) - emf_3d_in(i,j,k,3,3))             
                dxR_14 = idx2 * (emf_3d_in(i,j+1,k,3,1) - emf_3d_in(i,j+1,k,3,3))            
                dxL_34 = idx2 * (emf_3d_in(i+1,j,k,3,3) - emf_3d_in(i,j,k,3,1))
                dxR_34 = idx2 * (emf_3d_in(i+1,j+1,k,3,3) - emf_3d_in(i,j+1,k,3,1))
                        
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
             
                ! cell-corner electric field z component: Ez_i+1/2,j+1/2,k 
                emf_corner_in(i,j,k,3) =  0.25d0 * (emf_3d_in(i,j,k,3,1) + emf_3d_in(i,j+1,k,3,1) + &
                                       emf_3d_in(i,j,k,3,2) + emf_3d_in(i+1,j,k,3,2)) + &
                                       0.125d0 * dx * ((dy_12_14 - dy_12_34) + (dx_14_12 - dx_34_12)) 
                        
            END DO
        END DO    
    END DO    


    ! Corner electric field y-component: Ey_i+1/2,j,k+1/2
    DO k = klow-nb+1, khi+nb-2
        DO j = jlow-nb+1, jhi+nb-2
            DO i = ilow-nb+1, ihi+nb-2 
            
                ! get normal velocity at x-interface (roe-average value)
                vx = 0.5d0 * (q_3d(i,j,k,2)/q_3d(i,j,k,1) + q_3d(i+1,j,k,2)/q_3d(i+1,j,k,1) )
            
                ! get normal velocity at z-interface (roe-average value)
                vz = 0.5d0 * (q_3d(i,j,k,4)/q_3d(i,j,k,1) + q_3d(i,j,k+1,4)/q_3d(i,j,k+1,1) )
    
                ! compute upwinded interface gradients of electric field at x-interface
                dzL_14 = idx2 * (emf_3d_in(i,j,k,2,2) - emf_3d_in(i,j,k,2,3)) 
                dzR_14 = idx2 * (emf_3d_in(i+1,j,k,2,2) - emf_3d_in(i+1,j,k,2,3)) 
                dzL_34 = idx2 * (emf_3d_in(i,j,k+1,2,3) - emf_3d_in(i,j,k,2,2)) 
                dzR_34 = idx2 * (emf_3d_in(i+1,j,k+1,2,3) - emf_3d_in(i+1,j,k,2,2))
                
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
                dxL_14 = idx2 * (emf_3d_in(i,j,k,2,1) - emf_3d_in(i,j,k,2,3))             
                dxR_14 = idx2 * (emf_3d_in(i,j,k+1,2,1) - emf_3d_in(i,j,k+1,2,3))            
                dxL_34 = idx2 * (emf_3d_in(i+1,j,k,2,3) - emf_3d_in(i,j,k,2,1))
                dxR_34 = idx2 * (emf_3d_in(i+1,j,k+1,2,3) - emf_3d_in(i,j,k+1,2,1))
                        
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
             
                ! cell-corner electric field y component: Ey_i+1/2,j,k+1/2 
                emf_corner_in(i,j,k,2) =  0.25d0 * (emf_3d_in(i,j,k,2,1) + emf_3d_in(i,j,k+1,2,1) + &
                                       emf_3d_in(i,j,k,2,2) + emf_3d_in(i+1,j,k,2,2)) + &
                                       0.125d0 * dx * ((dz_12_14 - dz_12_34) + (dx_14_12 - dx_34_12)) 
                        
            END DO
        END DO    
    END DO    

    ! Corner electric field x-component: Ex_i,j+1/2,k+1/2
    DO k = klow-nb+1, khi+nb-2
        DO j = jlow-nb+1, jhi+nb-2
            DO i = ilow-nb+1, ihi+nb-2  
        
                ! get normal velocity at y-interface (roe-average value)
                vy = 0.5d0 * (q_3d(i,j,k,3)/q_3d(i,j,k,1) + q_3d(i,j+1,k,3)/q_3d(i,j+1,k,1) )
                
                ! get normal velocity at z-interface (roe-average value)
                vz = 0.5d0 * (q_3d(i,j,k,4)/q_3d(i,j,k,1) + q_3d(i,j,k+1,4)/q_3d(i,j,k+1,1) )
    
                ! compute upwinded interface gradients of electric field at y-interface
                dzL_14 = idx2 * (emf_3d_in(i,j,k,1,2) - emf_3d_in(i,j,k,1,3)) 
                dzR_14 = idx2 * (emf_3d_in(i,j+1,k,1,2) - emf_3d_in(i,j+1,k,1,3)) 
                dzL_34 = idx2 * (emf_3d_in(i,j,k+1,1,3) - emf_3d_in(i,j,k,1,2)) 
                dzR_34 = idx2 * (emf_3d_in(i,j+1,k+1,1,3) - emf_3d_in(i,j+1,k,1,2))
                
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
                dyL_14 = idx2 * (emf_3d_in(i,j,k,1,1) - emf_3d_in(i,j,k,1,3))             
                dyR_14 = idx2 * (emf_3d_in(i,j,k+1,1,1) - emf_3d_in(i,j,k+1,1,3))            
                dyL_34 = idx2 * (emf_3d_in(i,j+1,k,1,3) - emf_3d_in(i,j,k,1,1))
                dyR_34 = idx2 * (emf_3d_in(i,j+1,k+1,1,3) - emf_3d_in(i,j,k+1,1,1))
                        
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
             
                ! cell-corner electric field x component: Ex_i,j+1/2,k+1/2 
                emf_corner_in(i,j,k,1) =  0.25d0 * (emf_3d_in(i,j,k,1,1) + emf_3d_in(i,j,k+1,1,1) + &
                                       emf_3d_in(i,j,k,1,2) + emf_3d_in(i,j+1,k,1,2)) + &
                                       0.125d0 * dx * ((dz_12_14 - dz_12_34) + (dy_14_12 - dy_34_12)) 
                        
            END DO
        END DO    
    END DO    


END SUBROUTINE compute_cell_corner_emf_3d



! corrector stage half-dt update of state vector along x-direction using y and z fluxes
SUBROUTINE halfupdate_x_3d(Flux)

    REAL*8, INTENT(INOUT) :: Flux(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,7,3)
    INTEGER :: i,j,k
    REAL*8 :: dtdx, dBx, dBy, dBz, vy, vz, smomx, smomy, smomz, sBy, sBz

    dtdx = 0.5d0 * dt / dx

    IF(flux_protection) THEN

        protection_flags = .FALSE.
   
        ! make a copy of the initial state
        DO k = 1-nb+2, nz+nb-2
            DO j = 1-nb+2, ny+nb-2
                DO i = 1-nb+2, nx+nb-2
                    q_3d_prot(i,j,k,:) = q_3d(i,j,k,:)
                END DO
            END DO
        END DO  
        
    END IF

    ! Evolve the state vector by a half-stime step using only y and z fluxes
    DO k = 1-nb+2, nz+nb-2
        DO j = 1-nb+2, ny+nb-2
            DO i = 1-nb+2, nx+nb-2
   
                ! set the multi-dimensinal source terms
                dBx   = bface_work(i,j,k,1) - bface_work(i-1,j,k,1)   ! using initial cell-face magnetic field
                dBy   = bface_work(i,j,k,2) - bface_work(i,j-1,k,2)   ! using initial cell-face magnetic field
                dBz   = bface_work(i,j,k,3) - bface_work(i,j,k-1,3)   ! using initial cell-face magnetic field
                
                smomx = q_3d(i,j,k,5) * dBx
                smomy = q_3d(i,j,k,6) * dBx
                smomz = q_3d(i,j,k,7) * dBx
                
                vy = q_3d(i,j,k,3) / q_3d(i,j,k,1)
                vz = q_3d(i,j,k,4) / q_3d(i,j,k,1)
                
                sBy = 0.d0
                IF((-dBz * dBx) .GT. ZERO) sBy = -vy * SIGN(1.D0, dBz) * MIN( ABS(dBz), ABS(dBx) ) 
                
                sBz = 0.d0
                IF((-dBy * dBx) .GT. ZERO) sBz = -vz * SIGN(1.D0, dBy) * MIN( ABS(dBy), ABS(dBx) ) 
                                
                ! evolve using the y and z-fluxes
                q_3d(i,j,k,1) = q_3d(i,j,k,1) - dtdx * ( Flux(i,j,k,1,2) - Flux(i,j-1,k,1,2) + &
                                Flux(i,j,k,1,3) - Flux(i,j,k-1,1,3)) 
                q_3d(i,j,k,2) = q_3d(i,j,k,2) - dtdx * ( Flux(i,j,k,2,2) - Flux(i,j-1,k,2,2) + &
                                Flux(i,j,k,2,3) - Flux(i,j,k-1,2,3)) + dtdx * smomx
                q_3d(i,j,k,3) = q_3d(i,j,k,3) - dtdx * ( Flux(i,j,k,3,2) - Flux(i,j-1,k,3,2) + &
                                Flux(i,j,k,3,3) - Flux(i,j,k-1,3,3)) + dtdx * smomy
                q_3d(i,j,k,4) = q_3d(i,j,k,4) - dtdx * ( Flux(i,j,k,4,2) - Flux(i,j-1,k,4,2) + &
                                Flux(i,j,k,4,3) - Flux(i,j,k-1,4,3)) + dtdx * smomz 
                q_3d(i,j,k,5) = q_3d(i,j,k,5) - dtdx * ( Flux(i,j,k,5,2) - Flux(i,j-1,k,5,2) + &
                                Flux(i,j,k,5,3) - Flux(i,j,k-1,5,3))                                                              
                                
                q_3d(i,j,k,6) = q_3d(i,j,k,6) - 0.5d0 * dtdx * ( (emf_corner(i,j,k,1) - emf_corner(i,j,k-1,1)) + &
                                (emf_corner(i,j-1,k,1) - emf_corner(i,j-1,k-1,1)) ) + dtdx * sBy
                                
                q_3d(i,j,k,7) = q_3d(i,j,k,7) + 0.5d0 * dtdx * ( (emf_corner(i,j,k,1) - emf_corner(i,j-1,k,1)) + &
                                (emf_corner(i,j,k-1,1) - emf_corner(i,j-1,k-1,1)) ) + dtdx * sBz 

               
            END DO
        END DO
    END DO

    ! check for negative densities    
    IF(flux_protection) THEN

        DO k = 1-nb+2, nz+nb-2
            DO j = 1-nb+2, ny+nb-2
                DO i = 1-nb+2, nx+nb-2
   
                    IF(q_3d(i,j,k,1) .LT. MIN_DENS) THEN
                    
                        protection_flags(i,j,k) = .TRUE.
                        
                        ! copy protection fluxes through boundary of bad cell into main flux array
                        Flux(i,j,k,:,2)   = Fs_buff_prot(i,j,k,:,2)
                        Flux(i,j-1,k,:,2) = Fs_buff_prot(i,j-1,k,:,2)
                        Flux(i,j,k,:,3)   = Fs_buff_prot(i,j,k,:,3)
                        Flux(i,j,k-1,:,3) = Fs_buff_prot(i,j,k-1,:,3)
                        
                        emf_corner(i,j,k,1)     =  emf_corner_prot(i,j,k,1)
                        emf_corner(i,j-1,k,1)   =  emf_corner_prot(i,j-1,k,1)
                        emf_corner(i,j,k-1,1)   =  emf_corner_prot(i,j,k-1,1)
                        emf_corner(i,j-1,k-1,1) =  emf_corner_prot(i,j-1,k-1,1)
                        
                    END IF
   
                END DO
            END DO
        END DO
        
        ! apply protection fluxes to bad cells and their immediate neighbors
        ! Note: for the neighbor cells, the protection flux only applied to the
        ! boundary that is touching the bad cell.
        DO k = 1-nb+2, nz+nb-2
            DO j = 1-nb+2, ny+nb-2
                DO i = 1-nb+2, nx+nb-2
       
                    IF(protection_flags(i,j,k) .OR. protection_flags(i,j-1,k) .OR. protection_flags(i,j+1,k) .OR. &
                       protection_flags(i,j,k-1) .OR. protection_flags(i,j,k+1)) THEN
        
                        ! set the multi-dimensinal source terms
                        dBx   = bface_work(i,j,k,1) - bface_work(i-1,j,k,1)   ! using initial cell-face magnetic field
                        dBy   = bface_work(i,j,k,2) - bface_work(i,j-1,k,2)   ! using initial cell-face magnetic field
                        dBz   = bface_work(i,j,k,3) - bface_work(i,j,k-1,3)   ! using initial cell-face magnetic field
                        
                        smomx = q_3d_prot(i,j,k,5) * dBx
                        smomy = q_3d_prot(i,j,k,6) * dBx
                        smomz = q_3d_prot(i,j,k,7) * dBx
                        
                        vy = q_3d_prot(i,j,k,3) / q_3d_prot(i,j,k,1)
                        vz = q_3d_prot(i,j,k,4) / q_3d_prot(i,j,k,1)
                        
                        sBy = 0.d0
                        IF((-dBz * dBx) .GT. ZERO) sBy = -vy * SIGN(1.D0, dBz) * MIN( ABS(dBz), ABS(dBx) ) 
                        
                        sBz = 0.d0
                        IF((-dBy * dBx) .GT. ZERO) sBz = -vz * SIGN(1.D0, dBy) * MIN( ABS(dBy), ABS(dBx) ) 
                                        
                        ! evolve using the y and z-fluxes
                        q_3d(i,j,k,1) = q_3d_prot(i,j,k,1) - dtdx * ( Flux(i,j,k,1,2) - Flux(i,j-1,k,1,2) + &
                                        Flux(i,j,k,1,3) - Flux(i,j,k-1,1,3)) 
                        q_3d(i,j,k,2) = q_3d_prot(i,j,k,2) - dtdx * ( Flux(i,j,k,2,2) - Flux(i,j-1,k,2,2) + &
                                        Flux(i,j,k,2,3) - Flux(i,j,k-1,2,3)) + dtdx * smomx
                        q_3d(i,j,k,3) = q_3d_prot(i,j,k,3) - dtdx * ( Flux(i,j,k,3,2) - Flux(i,j-1,k,3,2) + &
                                        Flux(i,j,k,3,3) - Flux(i,j,k-1,3,3)) + dtdx * smomy
                        q_3d(i,j,k,4) = q_3d_prot(i,j,k,4) - dtdx * ( Flux(i,j,k,4,2) - Flux(i,j-1,k,4,2) + &
                                        Flux(i,j,k,4,3) - Flux(i,j,k-1,4,3)) + dtdx * smomz 
                        q_3d(i,j,k,5) = q_3d_prot(i,j,k,5) - dtdx * ( Flux(i,j,k,5,2) - Flux(i,j-1,k,5,2) + &
                                        Flux(i,j,k,5,3) - Flux(i,j,k-1,5,3))                                                              
                                        
                        q_3d(i,j,k,6) = q_3d_prot(i,j,k,6) - 0.5d0 * dtdx * ( (emf_corner(i,j,k,1) - emf_corner(i,j,k-1,1)) + &
                                        (emf_corner(i,j-1,k,1) - emf_corner(i,j-1,k-1,1)) ) + dtdx * sBy
                                        
                        q_3d(i,j,k,7) = q_3d_prot(i,j,k,7) + 0.5d0 * dtdx * ( (emf_corner(i,j,k,1) - emf_corner(i,j-1,k,1)) + &
                                        (emf_corner(i,j,k-1,1) - emf_corner(i,j-1,k-1,1)) ) + dtdx * sBz 

    
                    END IF
        
                END DO
            END DO
        END DO
        
    END IF
    
    CALL final_check_unphysical(1-nb+2, nx+nb-2, 1-nb+2, ny+nb-2, 1-nb+2, nz+nb-2)
    
END SUBROUTINE halfupdate_x_3d


! corrector stage half-dt update of state vector along y-direction using x-fluxes
SUBROUTINE halfupdate_y_3d(Flux)

    REAL*8, INTENT(INOUT) :: Flux(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,7,3)
    INTEGER :: i,j,k
    REAL*8 :: dtdx, dBx, dBy, dBz, vx, vz, smomx, smomy, smomz, sBx, sBz

    dtdx = 0.5d0 * dt / dx

    IF(flux_protection) THEN

        protection_flags = .FALSE.
   
        ! make a copy of the initial state
        DO k = 1-nb+2, nz+nb-2
            DO j = 1-nb+2, ny+nb-2
                DO i = 1-nb+2, nx+nb-2
                    q_3d_prot(i,j,k,:) = q_3d(i,j,k,:)
                END DO
            END DO
        END DO  
    
    END IF
    
    ! Evolve the state vector by a half-stime step using only y and z fluxes
    DO k = 1-nb+2, nz+nb-2
        DO j = 1-nb+2, ny+nb-2
            DO i = 1-nb+2, nx+nb-2
     
                ! set the multi-dimensinal source terms
                dBx   = bface_work(i,j,k,1) - bface_work(i-1,j,k,1)   ! using initial cell-face magnetic field
                dBy   = bface_work(i,j,k,2) - bface_work(i,j-1,k,2)   ! using initial cell-face magnetic field
                dBz   = bface_work(i,j,k,3) - bface_work(i,j,k-1,3)   ! using initial cell-face magnetic field
                
                smomx = q_3d(i,j,k,5) * dBy
                smomy = q_3d(i,j,k,6) * dBy
                smomz = q_3d(i,j,k,7) * dBy
                
                vx = q_3d(i,j,k,2) / q_3d(i,j,k,1)
                vz = q_3d(i,j,k,4) / q_3d(i,j,k,1)
                
                sBz = 0.d0
                IF((-dBx * dBy) .GT. ZERO) sBz = -vz * SIGN(1.D0, dBx) * MIN( ABS(dBx), ABS(dBy) ) 
                
                sBx = 0.d0
                IF((-dBz * dBy) .GT. ZERO) sBx = -vx * SIGN(1.D0, dBz) * MIN( ABS(dBz), ABS(dBy) ) 
                                
                ! evolve using the x and z-fluxes
                q_3d(i,j,k,1) = q_3d(i,j,k,1) - dtdx * ( Flux(i,j,k,1,1) - Flux(i-1,j,k,1,1) + &
                                Flux(i,j,k,1,3) - Flux(i,j,k-1,1,3)) 
                q_3d(i,j,k,2) = q_3d(i,j,k,2) - dtdx * ( Flux(i,j,k,2,1) - Flux(i-1,j,k,2,1) + &
                                Flux(i,j,k,2,3) - Flux(i,j,k-1,2,3)) + dtdx * smomx
                q_3d(i,j,k,3) = q_3d(i,j,k,3) - dtdx * ( Flux(i,j,k,3,1) - Flux(i-1,j,k,3,1) + &
                                Flux(i,j,k,3,3) - Flux(i,j,k-1,3,3)) + dtdx * smomy
                q_3d(i,j,k,4) = q_3d(i,j,k,4) - dtdx * ( Flux(i,j,k,4,1) - Flux(i-1,j,k,4,1) + &
                                Flux(i,j,k,4,3) - Flux(i,j,k-1,4,3)) + dtdx * smomz 
                                
                q_3d(i,j,k,5) = q_3d(i,j,k,5) + 0.5d0 * dtdx * ( (emf_corner(i,j,k,2) - emf_corner(i,j,k-1,2)) + &
                                (emf_corner(i-1,j,k,2) - emf_corner(i-1,j,k-1,2)) ) + dtdx * sBx
                                
                q_3d(i,j,k,6) = q_3d(i,j,k,6) - dtdx * ( Flux(i,j,k,6,1) - Flux(i-1,j,k,6,1) + &
                                Flux(i,j,k,6,3) - Flux(i,j,k-1,6,3))         
                                
                q_3d(i,j,k,7) = q_3d(i,j,k,7) - 0.5d0 * dtdx * ( (emf_corner(i,j,k,2) - emf_corner(i-1,j,k,2)) + &
                                (emf_corner(i,j,k-1,2) - emf_corner(i-1,j,k-1,2)) ) + dtdx * sBz


            END DO
        END DO
    END DO

    ! check for negative densities    
    IF(flux_protection) THEN

        DO k = 1-nb+2, nz+nb-2
            DO j = 1-nb+2, ny+nb-2
                DO i = 1-nb+2, nx+nb-2
   
                    IF(q_3d(i,j,k,1) .LT. MIN_DENS) THEN
                    
                        protection_flags(i,j,k) = .TRUE.
                        
                        ! copy protection fluxes through boundary of bad cell into main flux array
                        Flux(i,j,k,:,1)   = Fs_buff_prot(i,j,k,:,1)
                        Flux(i-1,j,k,:,1) = Fs_buff_prot(i-1,j,k,:,1)
                        Flux(i,j,k,:,3)   = Fs_buff_prot(i,j,k,:,3)
                        Flux(i,j,k-1,:,3) = Fs_buff_prot(i,j,k-1,:,3)
                        
                        emf_corner(i,j,k,2)     =  emf_corner_prot(i,j,k,2)
                        emf_corner(i-1,j,k,2)   =  emf_corner_prot(i-1,j,k,2)
                        emf_corner(i,j,k-1,2)   =  emf_corner_prot(i,j,k-1,2)
                        emf_corner(i-1,j,k-1,2) =  emf_corner_prot(i-1,j,k-1,2)

                    END IF
       
                END DO
            END DO
        END DO

        DO k = 1-nb+2, nz+nb-2
            DO j = 1-nb+2, ny+nb-2
                DO i = 1-nb+2, nx+nb-2
       
                    IF(protection_flags(i,j,k) .OR. protection_flags(i-1,j,k) .OR. protection_flags(i+1,j,k) .OR. &
                       protection_flags(i,j,k-1) .OR. protection_flags(i,j,k+1)) THEN
        
                        ! set the multi-dimensinal source terms
                        dBx   = bface_work(i,j,k,1) - bface_work(i-1,j,k,1)   ! using initial cell-face magnetic field
                        dBy   = bface_work(i,j,k,2) - bface_work(i,j-1,k,2)   ! using initial cell-face magnetic field
                        dBz   = bface_work(i,j,k,3) - bface_work(i,j,k-1,3)   ! using initial cell-face magnetic field
                        
                        smomx = q_3d_prot(i,j,k,5) * dBy
                        smomy = q_3d_prot(i,j,k,6) * dBy
                        smomz = q_3d_prot(i,j,k,7) * dBy
                        
                        vx = q_3d_prot(i,j,k,2) / q_3d_prot(i,j,k,1)
                        vz = q_3d_prot(i,j,k,4) / q_3d_prot(i,j,k,1)
                        
                        sBz = 0.d0
                        IF((-dBx * dBy) .GT. ZERO) sBz = -vz * SIGN(1.D0, dBx) * MIN( ABS(dBx), ABS(dBy) ) 
                        
                        sBx = 0.d0
                        IF((-dBz * dBy) .GT. ZERO) sBx = -vx * SIGN(1.D0, dBz) * MIN( ABS(dBz), ABS(dBy) ) 
                                        
                        ! evolve using the x and z-fluxes
                        q_3d(i,j,k,1) = q_3d_prot(i,j,k,1) - dtdx * ( Flux(i,j,k,1,1) - Flux(i-1,j,k,1,1) + &
                                        Flux(i,j,k,1,3) - Flux(i,j,k-1,1,3)) 
                        q_3d(i,j,k,2) = q_3d_prot(i,j,k,2) - dtdx * ( Flux(i,j,k,2,1) - Flux(i-1,j,k,2,1) + &
                                        Flux(i,j,k,2,3) - Flux(i,j,k-1,2,3)) + dtdx * smomx
                        q_3d(i,j,k,3) = q_3d_prot(i,j,k,3) - dtdx * ( Flux(i,j,k,3,1) - Flux(i-1,j,k,3,1) + &
                                        Flux(i,j,k,3,3) - Flux(i,j,k-1,3,3)) + dtdx * smomy
                        q_3d(i,j,k,4) = q_3d_prot(i,j,k,4) - dtdx * ( Flux(i,j,k,4,1) - Flux(i-1,j,k,4,1) + &
                                        Flux(i,j,k,4,3) - Flux(i,j,k-1,4,3)) + dtdx * smomz 
                                        
                        q_3d(i,j,k,5) = q_3d_prot(i,j,k,5) + 0.5d0 * dtdx * ( (emf_corner(i,j,k,2) - emf_corner(i,j,k-1,2)) + &
                                        (emf_corner(i-1,j,k,2) - emf_corner(i-1,j,k-1,2)) ) + dtdx * sBx
                                        
                        q_3d(i,j,k,6) = q_3d_prot(i,j,k,6) - dtdx * ( Flux(i,j,k,6,1) - Flux(i-1,j,k,6,1) + &
                                        Flux(i,j,k,6,3) - Flux(i,j,k-1,6,3))         
                                        
                        q_3d(i,j,k,7) = q_3d_prot(i,j,k,7) - 0.5d0 * dtdx * ( (emf_corner(i,j,k,2) - emf_corner(i-1,j,k,2)) + &
                                        (emf_corner(i,j,k-1,2) - emf_corner(i-1,j,k-1,2)) ) + dtdx * sBz


                    END IF
       
                END DO
            END DO
        END DO

    END IF

    CALL final_check_unphysical(1-nb+2, nx+nb-2, 1-nb+2, ny+nb-2, 1-nb+2, nz+nb-2)
    
    
END SUBROUTINE halfupdate_y_3d


! corrector stage half-dt update of state vector along y-direction using x-fluxes
SUBROUTINE halfupdate_z_3d(Flux)

    REAL*8, INTENT(INOUT) :: Flux(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,7,3)
    INTEGER :: i,j,k
    REAL*8 :: dtdx, dBx, dBy, dBz, vx, vy, smomx, smomy, smomz, sBx, sBy

    dtdx = 0.5d0 * dt / dx

    IF(flux_protection) THEN

        protection_flags = .FALSE.
   
        ! make a copy of the initial state
        DO k = 1-nb+2, nz+nb-2
            DO j = 1-nb+2, ny+nb-2
                DO i = 1-nb+2, nx+nb-2
                    q_3d_prot(i,j,k,:) = q_3d(i,j,k,:)
                END DO
            END DO
        END DO  
        
    END IF

    ! Evolve the state vector by a half-stime step using only y and z fluxes
    DO k = 1-nb+2, nz+nb-2
        DO j = 1-nb+2, ny+nb-2
            DO i = 1-nb+2, nx+nb-2
           
                ! set the multi-dimensinal source terms
                dBx   = bface_work(i,j,k,1) - bface_work(i-1,j,k,1)   ! using initial cell-face magnetic field
                dBy   = bface_work(i,j,k,2) - bface_work(i,j-1,k,2)   ! using initial cell-face magnetic field
                dBz   = bface_work(i,j,k,3) - bface_work(i,j,k-1,3)   ! using initial cell-face magnetic field
                
                smomx = q_3d(i,j,k,5) * dBz
                smomy = q_3d(i,j,k,6) * dBz
                smomz = q_3d(i,j,k,7) * dBz
                
                vx = q_3d(i,j,k,2) / q_3d(i,j,k,1)
                vy = q_3d(i,j,k,3) / q_3d(i,j,k,1)
                
                sBx = 0.d0
                IF((-dBy * dBz) .GT. ZERO) sBx = -vx * SIGN(1.D0, dBy) * MIN( ABS(dBy), ABS(dBz) ) 
                
                sBy = 0.d0
                IF((-dBx * dBz) .GT. ZERO) sBy = -vy * SIGN(1.D0, dBx) * MIN( ABS(dBx), ABS(dBz) ) 
                
                                
                ! evolve using the x and z-fluxes
                q_3d(i,j,k,1) = q_3d(i,j,k,1) - dtdx * ( Flux(i,j,k,1,1) - Flux(i-1,j,k,1,1) + &
                                Flux(i,j,k,1,2) - Flux(i,j-1,k,1,2)) 
                q_3d(i,j,k,2) = q_3d(i,j,k,2) - dtdx * ( Flux(i,j,k,2,1) - Flux(i-1,j,k,2,1) + &
                                Flux(i,j,k,2,2) - Flux(i,j-1,k,2,2)) + dtdx * smomx
                q_3d(i,j,k,3) = q_3d(i,j,k,3) - dtdx * ( Flux(i,j,k,3,1) - Flux(i-1,j,k,3,1) + &
                                Flux(i,j,k,3,2) - Flux(i,j-1,k,3,2)) + dtdx * smomy
                q_3d(i,j,k,4) = q_3d(i,j,k,4) - dtdx * ( Flux(i,j,k,4,1) - Flux(i-1,j,k,4,1) + &
                                Flux(i,j,k,4,2) - Flux(i,j-1,k,4,2)) + dtdx * smomz 
                                
                q_3d(i,j,k,5) = q_3d(i,j,k,5) - 0.5d0 * dtdx * ( (emf_corner(i,j,k,3) - emf_corner(i,j-1,k,3)) + &
                                (emf_corner(i-1,j,k,3) - emf_corner(i-1,j-1,k,3)) ) + dtdx * sBx
                                
                q_3d(i,j,k,6) = q_3d(i,j,k,6) + 0.5d0 * dtdx * ( (emf_corner(i,j,k,3) - emf_corner(i-1,j,k,3)) + &
                                (emf_corner(i,j-1,k,3) - emf_corner(i-1,j-1,k,3)) ) + dtdx * sBy               
                                
                q_3d(i,j,k,7) = q_3d(i,j,k,7)- dtdx * ( Flux(i,j,k,7,1) - Flux(i-1,j,k,7,1) + &
                                Flux(i,j,k,7,2) - Flux(i,j-1,k,7,2))
                                
            END DO
        END DO
    END DO
    
    
    ! check for unphysical state
    IF(flux_protection) THEN

        DO k = 1-nb+2, nz+nb-2
            DO j = 1-nb+2, ny+nb-2
                DO i = 1-nb+2, nx+nb-2
       
                    IF(q_3d(i,j,k,1) .LT. MIN_DENS) THEN
                    
                        protection_flags(i,j,k) = .TRUE.
                        
                        ! copy protection fluxes through boundary of bad cell into main flux array
                        Flux(i,j,k,:,1)   = Fs_buff_prot(i,j,k,:,1)
                        Flux(i-1,j,k,:,1) = Fs_buff_prot(i-1,j,k,:,1)
                        Flux(i,j,k,:,2)   = Fs_buff_prot(i,j,k,:,2)
                        Flux(i,j-1,k,:,2) = Fs_buff_prot(i,j-1,k,:,2)
                        
                        emf_corner(i,j,k,3)     =  emf_corner_prot(i,j,k,3)
                        emf_corner(i-1,j,k,3)   =  emf_corner_prot(i-1,j,k,3)
                        emf_corner(i,j-1,k,3)   =  emf_corner_prot(i,j-1,k,3)
                        emf_corner(i-1,j-1,k,3) =  emf_corner_prot(i-1,j-1,k,3)

                    END IF
       
                END DO
            END DO
        END DO
    
        DO k = 1-nb+2, nz+nb-2
            DO j = 1-nb+2, ny+nb-2
                DO i = 1-nb+2, nx+nb-2
   
                    IF(protection_flags(i,j,k) .OR. protection_flags(i-1,j,k) .OR. protection_flags(i+1,j,k) .OR. &
                       protection_flags(i,j-1,k) .OR. protection_flags(i,j+1,k)) THEN
                       
                        ! set the multi-dimensinal source terms
                        dBx   = bface_work(i,j,k,1) - bface_work(i-1,j,k,1)   ! using initial cell-face magnetic field
                        dBy   = bface_work(i,j,k,2) - bface_work(i,j-1,k,2)   ! using initial cell-face magnetic field
                        dBz   = bface_work(i,j,k,3) - bface_work(i,j,k-1,3)   ! using initial cell-face magnetic field
                        
                        smomx = q_3d_prot(i,j,k,5) * dBz
                        smomy = q_3d_prot(i,j,k,6) * dBz
                        smomz = q_3d_prot(i,j,k,7) * dBz
                        
                        vx = q_3d_prot(i,j,k,2) / q_3d_prot(i,j,k,1)
                        vy = q_3d_prot(i,j,k,3) / q_3d_prot(i,j,k,1)
                        
                        sBx = 0.d0
                        IF((-dBy * dBz) .GT. ZERO) sBx = -vx * SIGN(1.D0, dBy) * MIN( ABS(dBy), ABS(dBz) ) 
                        
                        sBy = 0.d0
                        IF((-dBx * dBz) .GT. ZERO) sBy = -vy * SIGN(1.D0, dBx) * MIN( ABS(dBx), ABS(dBz) ) 
                        
                                        
                        ! evolve using the x and z-fluxes
                        q_3d(i,j,k,1) = q_3d_prot(i,j,k,1) - dtdx * ( Flux(i,j,k,1,1) - Flux(i-1,j,k,1,1) + &
                                        Flux(i,j,k,1,2) - Flux(i,j-1,k,1,2)) 
                        q_3d(i,j,k,2) = q_3d_prot(i,j,k,2) - dtdx * ( Flux(i,j,k,2,1) - Flux(i-1,j,k,2,1) + &
                                        Flux(i,j,k,2,2) - Flux(i,j-1,k,2,2)) + dtdx * smomx
                        q_3d(i,j,k,3) = q_3d_prot(i,j,k,3) - dtdx * ( Flux(i,j,k,3,1) - Flux(i-1,j,k,3,1) + &
                                        Flux(i,j,k,3,2) - Flux(i,j-1,k,3,2)) + dtdx * smomy
                        q_3d(i,j,k,4) = q_3d_prot(i,j,k,4) - dtdx * ( Flux(i,j,k,4,1) - Flux(i-1,j,k,4,1) + &
                                        Flux(i,j,k,4,2) - Flux(i,j-1,k,4,2)) + dtdx * smomz 
                                        
                        q_3d(i,j,k,5) = q_3d_prot(i,j,k,5) - 0.5d0 * dtdx * ( (emf_corner(i,j,k,3) - emf_corner(i,j-1,k,3)) + &
                                        (emf_corner(i-1,j,k,3) - emf_corner(i-1,j-1,k,3)) ) + dtdx * sBx
                                        
                        q_3d(i,j,k,6) = q_3d_prot(i,j,k,6) + 0.5d0 * dtdx * ( (emf_corner(i,j,k,3) - emf_corner(i-1,j,k,3)) + &
                                        (emf_corner(i,j-1,k,3) - emf_corner(i-1,j-1,k,3)) ) + dtdx * sBy               
                                        
                        q_3d(i,j,k,7) = q_3d_prot(i,j,k,7)- dtdx * ( Flux(i,j,k,7,1) - Flux(i-1,j,k,7,1) + &
                                        Flux(i,j,k,7,2) - Flux(i,j-1,k,7,2))

                    END IF
       
                END DO
            END DO
        END DO
       
    END IF

    CALL final_check_unphysical(1-nb+2, nx+nb-2, 1-nb+2, ny+nb-2, 1-nb+2, nz+nb-2)


END SUBROUTINE halfupdate_z_3d



! final (unsplit) finite volume update of the state vector (excluding Bx and By, which require CT update)
SUBROUTINE update_state_vector_3d(Flux)

    REAL*8, INTENT(INOUT) :: Flux(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,7,3)
    REAL*8 :: dtdx
    INTEGER :: i, j, k
    
    dtdx = dt / dx
  
    IF(flux_protection) THEN

        protection_flags = .FALSE.
   
        ! make a copy of the initial state
        DO k = 1, nz     
            DO j = 1, ny     
                DO i = 1, nx
                    q_3d_prot(i,j,k,:) = q_3d(i,j,k,:)
                END DO
            END DO
        END DO     
    
    END IF
    
    ! update the conserved variables (density and momentum only)
    DO k = 1, nz     
        DO j = 1, ny     
            DO i = 1, nx
        
                q_3d(i,j,k,1) = q_3d(i,j,k,1) - dtdx * ( (Flux(i,j,k,1,1) - Flux(i-1,j,k,1,1)) + (Flux(i,j,k,1,2) - Flux(i,j-1,k,1,2)) + &
                                (Flux(i,j,k,1,3) - Flux(i,j,k-1,1,3)) ) 
                                
                q_3d(i,j,k,2) = q_3d(i,j,k,2) - dtdx * ( (Flux(i,j,k,2,1) - Flux(i-1,j,k,2,1)) + (Flux(i,j,k,2,2) - Flux(i,j-1,k,2,2)) + &
                                (Flux(i,j,k,2,3) - Flux(i,j,k-1,2,3)) )
                                
                q_3d(i,j,k,3) = q_3d(i,j,k,3) - dtdx * ( (Flux(i,j,k,3,1) - Flux(i-1,j,k,3,1)) + (Flux(i,j,k,3,2) - Flux(i,j-1,k,3,2)) + &
                                (Flux(i,j,k,3,3) - Flux(i,j,k-1,3,3)) )
                                
                q_3d(i,j,k,4) = q_3d(i,j,k,4) - dtdx * ( (Flux(i,j,k,4,1) - Flux(i-1,j,k,4,1)) + (Flux(i,j,k,4,2) - Flux(i,j-1,k,4,2)) + &
                                (Flux(i,j,k,4,3) - Flux(i,j,k-1,4,3)) )
                                
                !q_3d(i,j,k,5) = q_3d(i,j,k,5) - dtdx * ( (Flux(i,j,k,5,1) - Flux(i-1,j,k,5,1)) + (Flux(i,j,k,5,2) - Flux(i,j-1,k,5,2)) + &
                !                (Flux(i,j,k,5,3) - Flux(i,j,k-1,5,3)) )  
                                
                !q_3d(i,j,k,6) = q_3d(i,j,k,6) - dtdx * ( (Flux(i,j,k,6,1) - Flux(i-1,j,k,6,1)) + (Flux(i,j,k,6,2) - Flux(i,j-1,k,6,2)) + &
                !                (Flux(i,j,k,6,3) - Flux(i,j,k-1,6,3)) )  
                                
                !q_3d(i,j,k,7) = q_3d(i,j,k,7) - dtdx * ( (Flux(i,j,k,7,1) - Flux(i-1,j,k,7,1)) + (Flux(i,j,k,7,2) - Flux(i,j-1,k,7,2)) + &
                !                (Flux(i,j,k,7,3) - Flux(i,j,k-1,7,3)) )
                
            END DO    
        END DO
    END DO

    ! check for unphysical state
    IF(flux_protection) THEN

        DO k = 1, nz
            DO j = 1, ny
                DO i = 1, nx
       
                    IF(q_3d(i,j,k,1) .LT. MIN_DENS) THEN
                    
                        protection_flags(i,j,k) = .TRUE.
                        
                        ! copy protection fluxes through boundary of bad cell into main Flux array
                        Flux(i,j,k,:,1)   = Fs_prot(i,j,k,:,1)
                        Flux(i-1,j,k,:,1) = Fs_prot(i-1,j,k,:,1)
                        Flux(i,j,k,:,2)   = Fs_prot(i,j,k,:,2)
                        Flux(i,j-1,k,:,2) = Fs_prot(i,j-1,k,:,2)
                        Flux(i,j,k,:,3)   = Fs_prot(i,j,k,:,3)
                        Flux(i,j,k-1,:,3) = Fs_prot(i,j,k-1,:,3)
                        
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

        DO k = 1, nz
            DO j = 1, ny
                DO i = 1, nx

                     IF(protection_flags(i,j,k) .OR. protection_flags(i-1,j,k) .OR. protection_flags(i+1,j,k) .OR. &
                        protection_flags(i,j-1,k) .OR. protection_flags(i,j+1,k) .OR. protection_flags(i,j,k-1) .OR. protection_flags(i,j,k+1)) THEN      

                        q_3d(i,j,k,1) = q_3d_prot(i,j,k,1) - dtdx * ( (Flux(i,j,k,1,1) - Flux(i-1,j,k,1,1)) + (Flux(i,j,k,1,2) - Flux(i,j-1,k,1,2)) + &
                                        (Flux(i,j,k,1,3) - Flux(i,j,k-1,1,3)) ) 
                                        
                        q_3d(i,j,k,2) = q_3d_prot(i,j,k,2) - dtdx * ( (Flux(i,j,k,2,1) - Flux(i-1,j,k,2,1)) + (Flux(i,j,k,2,2) - Flux(i,j-1,k,2,2)) + &
                                        (Flux(i,j,k,2,3) - Flux(i,j,k-1,2,3)) )
                                        
                        q_3d(i,j,k,3) = q_3d_prot(i,j,k,3) - dtdx * ( (Flux(i,j,k,3,1) - Flux(i-1,j,k,3,1)) + (Flux(i,j,k,3,2) - Flux(i,j-1,k,3,2)) + &
                                        (Flux(i,j,k,3,3) - Flux(i,j,k-1,3,3)) )
                                        
                        q_3d(i,j,k,4) = q_3d_prot(i,j,k,4) - dtdx * ( (Flux(i,j,k,4,1) - Flux(i-1,j,k,4,1)) + (Flux(i,j,k,4,2) - Flux(i,j-1,k,4,2)) + &
                                        (Flux(i,j,k,4,3) - Flux(i,j,k-1,4,3)) )
                                        
                        !q_3d(i,j,k,5) = q_3d(i,j,k,5) - dtdx * ( (Flux(i,j,k,5,1) - Flux(i-1,j,k,5,1)) + (Flux(i,j,k,5,2) - Flux(i,j-1,k,5,2)) + &
                        !                (Flux(i,j,k,5,3) - Flux(i,j,k-1,5,3)) )  
                                        
                        !q_3d(i,j,k,6) = q_3d(i,j,k,6) - dtdx * ( (Flux(i,j,k,6,1) - Flux(i-1,j,k,6,1)) + (Flux(i,j,k,6,2) - Flux(i,j-1,k,6,2)) + &
                        !                (Flux(i,j,k,6,3) - Flux(i,j,k-1,6,3)) )  
                                        
                        !q_3d(i,j,k,7) = q_3d(i,j,k,7) - dtdx * ( (Flux(i,j,k,7,1) - Flux(i-1,j,k,7,1)) + (Flux(i,j,k,7,2) - Flux(i,j-1,k,7,2)) + &
                        !                (Flux(i,j,k,7,3) - Flux(i,j,k-1,7,3)) )

                        
                    END IF
       
                END DO
            END DO
        END DO

    END IF

    CALL final_check_unphysical(1, nx, 1, ny, 1, nz)


END SUBROUTINE update_state_vector_3d


! updates cell-faces magnetic fields via CT
SUBROUTINE update_CT_3d(ilow, ihi, jlow, jhi, klow, khi, corrector_stage)

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi, klow, khi
    LOGICAL, INTENT(IN) :: corrector_stage
    REAL*8 :: dtdx
    INTEGER :: i, j, k

    dtdx = dt / dx

    ! for the initial predictor stage, we only apply a half-time update
    IF(.NOT. corrector_stage) dtdx = 0.5d0 * dtdx
    
    
    DO k = klow-nb+2, khi+nb-3
        DO j = jlow-nb+2, jhi+nb-3
            DO i = ilow-nb+2, ihi+nb-3  

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


    ! If this is the final (i.e. corrector stage) CT update, then also compute the cell-centered magnetic fields
    ! by averaging over the cell-face values (i.e. linear interpolation)
    IF(corrector_stage) THEN
    
        DO k = 1, nz
            DO j = 1, ny
                DO i = 1, nx  

                    ! Bx_i,j,k
                    q_3d(i,j,k,5) = 0.5d0 * (bface_3d(i,j,k,1) + bface_3d(i-1,j,k,1))  
                
                    ! By_i,j,k
                    q_3d(i,j,k,6) = 0.5d0 * (bface_3d(i,j,k,2) + bface_3d(i,j-1,k,2))

                    ! Bz_i,j,k
                    q_3d(i,j,k,7) = 0.5d0 * (bface_3d(i,j,k,3) + bface_3d(i,j,k-1,3))                    
               
                END DO
            END DO
        END DO
 
    END IF

END SUBROUTINE update_CT_3d



! Compute cell-interface roe averaged state
SUBROUTINE compute_roe_avg_state_1d(ilow, ihi, jlow, jhi)

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi
    REAL*8 :: rho2, vx2, vy2, vz2, Bx2, By2, Bz2, isrho
    INTEGER :: i, j
    
       
    DO j = jlow-nb, jhi+nb     
    DO i = ilow-nb, ihi+nb-1

        ! compute Roe-averaged primitive variables (this is just an arithmetic average)
        rho2 = 0.5d0 * ( q(i,j,1) + q(i+1,j,1) )
        isrho = 1.D0 / SQRT(rho2)
        vx2  = 0.5d0 * ( q(i,j,2)/q(i,j,1) + q(i+1,j,2)/q(i+1,j,1) )
        vy2  = 0.5d0 * ( q(i,j,3)/q(i,j,1) + q(i+1,j,3)/q(i+1,j,1) )
        vz2  = 0.5d0 * ( q(i,j,4)/q(i,j,1) + q(i+1,j,4)/q(i+1,j,1) )
        Bx2  = q(i,j,5) ! 0.5d0 * ( q(i,j,5) + q(i+1,j,5) )
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
        

    asqr  = sound_speed**2


    DO j = jlow-nb, jhi+nb     
    DO i = ilow-nb, ihi+nb-1

        rho2  = qavg(i,j,1)
        vx2   = qavg(i,j,2)
        bx    = qavg(i,j,5)
        bsqr  = qavg(i,j,5)**2 + qavg(i,j,6)**2 + qavg(i,j,7)**2
        btsqr  =  qavg(i,j,6)**2 + qavg(i,j,7)**2
        tmp1  = asqr + bsqr        
        tmp2  = SQRT( (asqr - bsqr)**2 + 4.d0 * asqr * btsqr )
        !tmp2  = SQRT(MAX(0.D0, tmp1**2 - 4.d0 * asqr * bx**2) )

        cs = SQRT( MAX(0.D0, 0.5d0 * (tmp1 - tmp2) ))  ! slow mode speed
        ca = ABS(bx)                                   ! alfven speed
        cf = SQRT(0.5d0 * (tmp1 + tmp2))               ! fast mode speed
        
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
    INTEGER :: i, j, k 
     

    asqr   = sound_speed**2      

         
    DO j = jlow-nb, jhi+nb     
    DO i = ilow-nb, ihi+nb-1

        bt = SQRT(qavg(i,j,6)**2 + qavg(i,j,7)**2) ! magnitude of B-transverse

        ! compute eigenvector renormalization co-efficients, taking limiting values for degenracies (K98 eqn. 2.13-15)
        
        !*************************** ORIGINAL *************************************
        !  Fails if c_s = c_a = c_f = c_sound
        !
        !IF(ABS(wave_speeds(i,j,3)**2 - wave_speeds(i,j,1)**2) .LT. ZERO) THEN
        !    alphaf(i,j) = 1.D0
        !    alphas(i,j) = 1.D0
        !ELSE
        !    tmp1 = 1.D0 / SQRT( MAX(0.d0, wave_speeds(i,j,3)**2 - wave_speeds(i,j,1)**2) )
        !    alphaf(i,j) = SQRT( MAX(0.d0, wave_speeds(i,j,3)**2 - wave_speeds(i,j,2)**2) ) * tmp1      ! Need to be careful with square roots. Sometimes round-off errors can make cf < ca, or ca < cs.   
        !    alphas(i,j) = SQRT( MAX(0.d0, wave_speeds(i,j,3)**2 - asqr) ) * tmp1
        !END IF
        
        
        !************************** MODIFIED V1 ********************************        
        !IF(ABS(wave_speeds(i,j,3)**2 - wave_speeds(i,j,1)**2) .LT. ZERO) THEN
        !    alphaf(i,j) = 1.D0
        !    alphas(i,j) = 1.D0
        !ELSE IF(ABS(wave_speeds(i,j,3)**2 - wave_speeds(i,j,2)**2) .LT. ZERO) THEN
        !    alphaf(i,j) = 0.D0
        !    alphas(i,j) = 1.D0
        !ELSE IF(ABS(wave_speeds(i,j,3)**2 - asqr) .LT. ZERO) THEN
        !    alphaf(i,j) = 1.D0
        !    alphas(i,j) = 0.D0            
        !ELSE
        !    tmp1 = 1.D0 / SQRT( MAX(0.d0, wave_speeds(i,j,3)**2 - wave_speeds(i,j,1)**2) )
        !    alphaf(i,j) = SQRT( MAX(0.d0, wave_speeds(i,j,3)**2 - wave_speeds(i,j,2)**2) ) * tmp1      ! Need to be careful with square roots. Sometimes round-off errors can make cf < ca, or ca < cs.   
        !    alphas(i,j) = SQRT( MAX(0.d0, wave_speeds(i,j,3)**2 - asqr) ) * tmp1
        !END IF
        !*************************************************************************
         
     
        !************************** MODIFIED V2 ********************************        
        IF((bt .LT. ZERO .AND. ABS(wave_speeds(i,j,2)**2 - asqr) .LT. ZERO) .OR. &
           (ABS(wave_speeds(i,j,3)**2 - wave_speeds(i,j,1)**2) .LT. ZERO) ) THEN        
            alphaf(i,j) = 1.D0
            alphas(i,j) = 1.D0   
        ELSE
            tmp1 = 1.D0 / SQRT( MAX(0.d0, wave_speeds(i,j,3)**2 - wave_speeds(i,j,1)**2) )
            alphaf(i,j) = SQRT( MAX(0.d0, wave_speeds(i,j,3)**2 - wave_speeds(i,j,2)**2) ) * tmp1      ! Need to be careful with square roots. Sometimes round-off errors can make cf < ca, or ca < cs.   
            alphas(i,j) = SQRT( MAX(0.d0, wave_speeds(i,j,3)**2 - asqr) ) * tmp1 
        END IF
        !*************************************************************************

        
        !IF(alphaf(i,j) .LT. ZERO .AND. alphas(i,j) .LT. ZERO) THEN
        !    PRINT*,'alphaf, alphas =',alphaf(i,j), alphas(i,j)
        !    PRINT*,'bt, ca^2-a^2, cf^2-cs^2 = ', bt, wave_speeds(i,j,2)**2 - asqr, wave_speeds(i,j,3)**2 - wave_speeds(i,j,1)**2
        !    STOP
        !END IF


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
    
    !*************************************************
    ! compute the right eigenvectors (K98 eq 2.12)
    !*************************************************
    
    DO j = jlow-nb, jhi+nb     
    DO i = ilow-nb, ihi+nb-1
    
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

    DO j = jlow-nb, jhi+nb     
    DO i = ilow-nb, ihi+nb-1
    
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
    
    DO j = jlow-nb, jhi+nb     
    DO i = ilow-nb, ihi+nb-1
    
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
        
    !*****************************************************
    ! compute the characteristic variables (K98 eqn. 2.27)
    !*****************************************************
    
    DO j = jlow-nb, jhi+nb     
    DO i = ilow-nb, ihi+nb-1
                            
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
    

    DO j = jlow-nb, jhi+nb     
    DO i = ilow-nb, ihi+nb-1
              
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
    
        
    DO j = jlow-nb, jhi+nb     
    DO i = ilow-nb, ihi+nb-1
              
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
        

END SUBROUTINE compute_eigenvectors_1d


! Compute the time-averaged cell-interface fluxes
SUBROUTINE compute_TVD_fluxes_1d(Flux, ilow, ihi, jlow, jhi, dtdx)

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi
    REAL*8, INTENT(IN) :: dtdx
    REAL*8, INTENT(INOUT) :: Flux(1-nb:nmax+nb,1-nb:nmax+nb,7)
    INTEGER :: i, j, k
    REAL*8 :: dxdt, gammak, betak
    REAL*8 :: Qk, chi
    REAL*8 :: gtilde1, gtilde2, Qk1, Qk2, chi1, chi2, signg
    REAL*8 :: gk(1-nb:nmax+nb, 1-nb:nmax+nb, nwaves), F(1-nb:ihi+nb, 1-nb:jhi+nb,nwaves), FTVD(1-nb:ihi+nb, 1-nb:jhi+nb,nwaves)
    REAL*8 :: rho, vx, vy, vz, Bx, By, Bz


    dxdt = dx/dt
    
    ! clear temp array
    FTVD = 0.D0
    
    ! compute cell center fluxes (instead of storing this in an array, maybe compute these ad hoc inside the TVD flux calculation loop?)
    DO j = jlow-nb, jhi+nb     
    DO i = ilow-nb, ihi+nb
    
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
        DO j = jlow-nb, jhi+nb     
        DO i = ilow-nb+1, ihi+nb-1
                
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
        DO j = jlow-nb, jhi+nb     
        DO i = ilow-nb+1, ihi+nb-2
        
        
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
            
            ! accumulate TVD FTVD contribution from each wave  (K98 eqn. 2.25)
            FTVD(i,j,1) = FTVD(i,j,1) - betak * Rk(i,j,1,k) 
            FTVD(i,j,2) = FTVD(i,j,2) - betak * Rk(i,j,2,k) 
            FTVD(i,j,3) = FTVD(i,j,3) - betak * Rk(i,j,3,k) 
            FTVD(i,j,4) = FTVD(i,j,4) - betak * Rk(i,j,4,k) 
            FTVD(i,j,5) = FTVD(i,j,5) - betak * Rk(i,j,5,k) 
            FTVD(i,j,6) = FTVD(i,j,6) - betak * Rk(i,j,6,k) 
            
        END DO    
        END DO
    END DO

    ! add remaining term in the Roe TVD flux
    DO j = jlow-nb, jhi+nb     
    DO i = ilow-nb+1, ihi+nb-2
        
        Flux(i,j,1) = 0.5d0 * ( dxdt * FTVD(i,j,1) + F(i,j,1) + F(i+1,j,1) )
        Flux(i,j,2) = 0.5d0 * ( dxdt * FTVD(i,j,2) + F(i,j,2) + F(i+1,j,2) )
        Flux(i,j,3) = 0.5d0 * ( dxdt * FTVD(i,j,3) + F(i,j,3) + F(i+1,j,3) )
        Flux(i,j,4) = 0.5d0 * ( dxdt * FTVD(i,j,4) + F(i,j,4) + F(i+1,j,4) )
        Flux(i,j,5) = 0.d0
        Flux(i,j,6) = 0.5d0 * ( dxdt * FTVD(i,j,5) + F(i,j,5) + F(i+1,j,5) )
        Flux(i,j,7) = 0.5d0 * ( dxdt * FTVD(i,j,6) + F(i,j,6) + F(i+1,j,6) )   
        
    END DO
    END DO
    
    
END SUBROUTINE compute_TVD_fluxes_1d



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
            uL(:) = q(i,j,:) 
            uR(:) = q(i+1,j,:) 
            uavg(:) = 0.5d0 * (uL(:) + uR(:)) ! Roe-averaged state
            
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
            tmp2 = SQRT( (asqr - bsqr)**2 + 4.d0 * asqr * btsqr )
            cfL   = SQRT( 0.5d0 * (tmp1 + tmp2) )
            
            lambdaf_minus_L = vxL - cfL  
            lambdaf_plus_L  = vxL + cfL  
            
            bsqr   = (uR(5)**2 + uR(6)**2 + uR(7)**2) / uR(1)
            btsqr  = (uR(6)**2 + uR(7)**2) / uR(1)
            
            tmp1 = asqr + bsqr        
            tmp2 = SQRT( (asqr - bsqr)**2 + 4.d0 * asqr * btsqr )
            cfR   = SQRT( 0.5d0 * (tmp1 + tmp2) )

            lambdaf_minus_R = vxR - cfR  
            lambdaf_plus_R  = vxR + cfR          


            bsqr  = (uavg(5)**2 + uavg(6)**2 + uavg(7)**2) / uavg(1)
            btsqr = (uavg(6)**2 + uavg(7)**2) / uavg(1)
            tmp1 = asqr + bsqr        
            tmp2 = SQRT( (asqr - bsqr)**2 + 4.d0 * asqr * btsqr )
            cf_avg   = SQRT( 0.5d0 * (tmp1 + tmp2) )
            
            vx_avg = uavg(2)/ uavg(1)

            lambdaf_minus_avg = vx_avg - cf_avg 
            lambdaf_plus_avg  = vx_avg + cf_avg

            sL = MIN(MIN(lambdaf_minus_L, lambdaf_minus_avg), 0.d0) !, 0.d0)
            sR = MAX(MAX(lambdaf_plus_R, lambdaf_plus_avg), 0.D0)
           
             
            isRL = 1.d0 / (sR - sL)

            ! compute HLL state
            uHLL(:) = (sR * uR(:) - sL * uL(:) - FR(:) + FL(:) ) * isRL
             
            ! compute HLL flux
            FHLL(:)  = (sR * FL(:) - sL * FR(:) + sR * sL* (uR(:) - uL(:)) ) * isRL


            Flux(i,j,:) = FHLL(:)               
           
            
        END DO    
    END DO   
    
    
    
END SUBROUTINE compute_HLLE_fluxes_1d



! this subroutine checks state vector array for negative densities 
SUBROUTINE final_check_unphysical(ilow, ihi, jlow, jhi, klow, khi)

    INTEGER, INTENT(IN) :: ilow, ihi, jlow, jhi, klow, khi
    INTEGER :: i, j, k

    DO k = klow, khi
        DO j = jlow, jhi
            DO i = ilow, ihi
    
                ! replace negative densities with user-defined floor value
                IF(q_3d(i,j,k,1) .LT. MIN_DENS) THEN
                    q_3d(i,j,k,1) = MIN_DENS 
                    !IF(print_debug) PRINT*,'Negative density detected in cell# ',i,j,k
                    PRINT*,'Negative density detected in cell# ',i,j,k
                END IF
                    
            END DO
        END DO
    END DO

END SUBROUTINE final_check_unphysical


! compute maximum speed (need this for calculating the time step size)
SUBROUTINE get_max_speed_tvd(max_speed)

    REAL*8, INTENT(INOUT) :: max_speed
    INTEGER :: i, j, k
    REAL*8 :: rho2, isrho, vx2, vy2, vz2, bx2, by2, bz2, bsqr, asqr, btsqr, cfx, cfy, cfz, tmp
    
    max_speed = 0.d0

    ! find the maximum fast mode speed on the grid (the fast mode speed at a cell-interface
    ! is computed using the roe-averaged state at that interface)
    DO k = 0, nz    
      DO j = 0, ny    
        DO i = 0, nx ! no need to include boundary cells for this

            asqr  = sound_speed**2
            
            rho2 = 0.5d0 * ( q_3d(i,j,k,1) + q_3d(i+1,j,k,1) )
            isrho = 1.D0 / SQRT(rho2)
            bx2  = 0.5d0 * ( q_3d(i,j,k,5) + q_3d(i+1,j,k,5) ) * isrho
            by2  = 0.5d0 * ( q_3d(i,j,k,6) + q_3d(i+1,j,k,6) ) * isrho
            bz2  = 0.5d0 * ( q_3d(i,j,k,7) + q_3d(i+1,j,k,7) ) * isrho
            btsqr = by2**2 + bz2**2 
            bsqr  = bx2**2 + btsqr 
            tmp = SQRT( (asqr - bsqr)**2 + 4.d0 * asqr * btsqr )
            
            cfx = SQRT( 0.5d0 * (asqr + bsqr + tmp) ) ! fast mode speed at x-interface

            rho2 = 0.5d0 * ( q_3d(i,j,k,1) + q_3d(i,j+1,k,1) )
            isrho = 1.D0 / SQRT(rho2)
            bx2  = 0.5d0 * ( q_3d(i,j,k,5) + q_3d(i,j+1,k,5) ) * isrho
            by2  = 0.5d0 * ( q_3d(i,j,k,6) + q_3d(i,j+1,k,6) ) * isrho
            bz2  = 0.5d0 * ( q_3d(i,j,k,7) + q_3d(i,j+1,k,7) ) * isrho
            btsqr = bx2**2 + bz2**2 
            bsqr  = by2**2 + btsqr 
            tmp = SQRT( (asqr - bsqr)**2 + 4.d0 * asqr * btsqr )
           
            cfy = SQRT( 0.5d0 * (asqr + bsqr + tmp) ) ! fast mode speed at y-interface

            rho2 = 0.5d0 * ( q_3d(i,j,k,1) + q_3d(i,j,k+1,1) )
            isrho = 1.D0 / SQRT(rho2)
            bx2  = 0.5d0 * ( q_3d(i,j,k,5) + q_3d(i,j,k+1,5) ) * isrho
            by2  = 0.5d0 * ( q_3d(i,j,k,6) + q_3d(i,j,k+1,6) ) * isrho
            bz2  = 0.5d0 * ( q_3d(i,j,k,7) + q_3d(i,j,k+1,7) ) * isrho
            btsqr = bx2**2 + by2**2 
            bsqr  = bz2**2 + btsqr 
            tmp = SQRT( (asqr - bsqr)**2 + 4.d0 * asqr * btsqr )
            
            cfz = SQRT( 0.5d0 * (asqr + bsqr + tmp) ) ! fast mode speed at z-interface
            
            vx2 = 0.5D0 *( (q_3d(i,j,k,2) /q_3d(i,j,k,1)) + (q_3d(i+1,j,k,2) /q_3d(i+1,j,k,1)) )
            vy2 = 0.5D0 *( (q_3d(i,j,k,3) /q_3d(i,j,k,1)) + (q_3d(i,j+1,k,3) /q_3d(i,j+1,k,1)) )
            vz2 = 0.5D0 *( (q_3d(i,j,k,4) /q_3d(i,j,k,1)) + (q_3d(i,j,k+1,4) /q_3d(i,j,k+1,1)) )

            max_speed = MAX(max_speed, ABS(vx2) + cfx, ABS(vy2) + cfy, ABS(vz2) + cfz)
            
        END DO    
      END DO
    END DO


    IF(max_speed .LT. ZERO) THEN
        PRINT*,'ERROR!! ZERO MAX SPEED!'
        STOP
    END IF


END SUBROUTINE get_max_speed_tvd


END MODULE Iso_TVD_Solver_mod