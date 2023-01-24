PROGRAM isothermal_tvd_driver

USE constants_mod
USE grid_data_mod
USE Iso_TVD_Solver_mod

IMPLICIT NONE

INTEGER :: tsteps


OPEN(UNIT=10, FILE='Output/fluid_xcut.txt')
OPEN(UNIT=11, FILE='Output/Bfield_xcut.txt')
OPEN(UNIT=12, FILE='Output/fluid_ycut.txt')
OPEN(UNIT=13, FILE='Output/Bfield_ycut.txt')



! allocate memory for grid data
CALL create_grid_arrays()

! initialize solver
CALL init_solver()

! initialize fluid state
CALL init_state()


! simulation loop
DO tsteps = 1, maxsteps

    PRINT*,'Time step = ',tsteps

    ! evolve state vector
    CALL solve_2d()
    
    ! apply boundary conditions
    CALL boundary_conditions()

    ! save state to file 
    IF(MOD(tsteps,tSkip) .EQ. 0) CALL fileOutput()

END DO


! ********************************
CALL destroy_solver()
CALL destroy_grid_arrays()
CLOSE(UNIT=10)
CLOSE(UNIT=11)

PRINT*,'Done!'

CONTAINS


SUBROUTINE init_state()

    INTEGER :: i, j
    REAL*8 :: rhoL, rhoR, vxL, vxR, vyL, vyR, vzL, vzR, &
              bxL, bxR, byL, byR, bzL, bzR
    
    
    ! Riemann problem left state
    rhoL = 1.d0
    vxL = 0.d0
    vyL = 0.d0
    vzL = 0.d0
    bxL = 0.d0/SQRT(FOURPI)
    byL = 5.d0/SQRT(FOURPI)
    bzL = 0.d0/SQRT(FOURPI)
    
    ! Riemann problem right state
    rhoR = 0.1d0
    vxR = 0.d0
    vyR = 0.d0
    vzR = 0.d0
    bxR = bxL
    byR = 2.d0/SQRT(FOURPI)
    bzR = 0.d0/SQRT(FOURPI)
    
    !GO TO 111
    
    DO i = 1-nb, nx+nb
                   
        IF(i .LT. nx/2) THEN
        
            q_2d(i,:,1) = rhoL
            q_2d(i,:,2) = rhoL*vxL        
            q_2d(i,:,3) = rhoL*vyL
            q_2d(i,:,4) = rhoL*vzL    
            q_2d(i,:,7) = bzL     
            
            bface_2d(i,:,1) = bxL
            bface_2d(i,:,2) = byL
            
        ELSE
            
            q_2d(i,:,1) = rhoR
            q_2d(i,:,2) = rhoR*vxR        
            q_2d(i,:,3) = rhoR*vyR
            q_2d(i,:,4) = rhoR*vzR    
            q_2d(i,:,7) = bzR     
           
            bface_2d(i,:,1) = bxR
            bface_2d(i,:,2) = byR
           
        END IF 
                
    END DO
    

    DO j = 1-nb+1, ny+nb    
    DO i = 1-nb+1, nx+nb

        q_2d(i,j,5) = 0.5d0 * ( bface_2d(i-1,j,1) + bface_2d(i,j,1) )     
        q_2d(i,j,6) = 0.5d0 * ( bface_2d(i,j-1,1) + bface_2d(i,j,1) )       
        
    END DO
    END DO
    
    !111 CONTINUE
    
    GO TO 112
    
    DO j = 1-nb, ny+nb
                   
        IF(j .LT. ny/2) THEN
        
            q_2d(:,j,1) = rhoL
            q_2d(:,j,2) = rhoL*vxL        
            q_2d(:,j,3) = rhoL*vyL
            q_2d(:,j,4) = rhoL*vzL    
            q_2d(:,j,5) = bxL    
            q_2d(:,j,6) = byL    
            q_2d(:,j,7) = bzL     
            
            bface_2d(:,j,1) = bxL
            bface_2d(:,j,2) = byL
            
        ELSE
            
            q_2d(:,j,1) = rhoR
            q_2d(:,j,2) = rhoR*vxR        
            q_2d(:,j,3) = rhoR*vyR
            q_2d(:,j,4) = rhoR*vzR    
            q_2d(:,j,5) = bxR    
            q_2d(:,j,6) = byR    
            q_2d(:,j,7) = bzR     
           
            bface_2d(:,j,1) = bxR
            bface_2d(:,j,2) = byR
           
        END IF 
                
    END DO
        
    112 CONTINUE
    
    
    CALL fileOutput()

    
END SUBROUTINE init_state



SUBROUTINE fileOutput()

    INTEGER :: i, j
    REAL*8 :: dx, vx, vy, vz
    
    dx = 1.d0/nx
    j = ny/2
    
    DO i =1,nx
        vx = q_2d(i,j,2)/q_2d(i,j,1)
        vy = q_2d(i,j,3)/q_2d(i,j,1)
        vz = q_2d(i,j,4)/q_2d(i,j,1)
        WRITE(10,*) i*dx,q_2d(i,j,1),vx,vy,vz
        WRITE(11,*) i*dx,q_2d(i,j,6),q_2d(i,j,7)        
    END DO

    i = nx/2
    
    DO j =1,ny
        vx = q_2d(i,j,2)/q_2d(i,j,1)
        vy = q_2d(i,j,3)/q_2d(i,j,1)
        vz = q_2d(i,j,4)/q_2d(i,j,1)
        WRITE(12,*) j*dx,q_2d(i,j,1),vx,vy,vz
        WRITE(13,*) j*dx,q_2d(i,j,6),q_2d(i,j,7)        
    END DO


END SUBROUTINE fileOutput


SUBROUTINE boundary_conditions()

     INTEGER :: i, j
     
    ! Open Boundaries
    
    DO j=1-nb,ny+nb
        q_2d(-3,j,:) = q_2d(-2,j,:)       
        q_2d(nx+2,j,:) = q_2d(nx+1,j,:)
        bface_2d(nx+1,j,:) = bface_2d(nx,j,:)  ! Is this correct?
        bface_2d(nx+2,j,:) = bface_2d(nx,j,:)  ! Is this correct?
        bface_2d(-1,j,:) = bface_2d(0,j,:)     ! ?
        bface_2d(-2,j,:) = bface_2d(0,j,:)     ! ?
    END DO

    DO i=1-nb,nx+nb
        q_2d(i,-3,:) = q_2d(i,-2,:)       
        q_2d(i,ny+2,:) = q_2d(i,ny+1,:)
        bface_2d(i,ny+1,:) = bface_2d(i,ny,:) ! ?
        bface_2d(i,ny+2,:) = bface_2d(i,ny,:) ! ?
        bface_2d(i,-1,:) = bface_2d(i,0,:)    ! ?
        bface_2d(i,-2,:) = bface_2d(i,0,:)    ! ?
    END DO


END SUBROUTINE boundary_conditions


END PROGRAM isothermal_tvd_driver