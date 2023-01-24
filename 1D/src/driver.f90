PROGRAM isothermal_tvd_driver

USE constants_mod
USE grid_data_mod
USE Iso_TVD_Solver_mod

IMPLICIT NONE

INTEGER :: tsteps


OPEN(UNIT=10, FILE='Output/fluid.txt')
OPEN(UNIT=11, FILE='Output/Bfield.txt')



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
    CALL solve()
    
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

    INTEGER :: i
    REAL*8 :: rhoL, rhoR, vxL, vxR, vyL, vyR, vzL, vzR, &
              bxL, bxR, byL, byR, bzL, bzR
    
    
    ! Riemann problem left state
    rhoL = 1.d0
    vxL = -1.d0
    vyL = 0.d0
    vzL = 0.d0
    bxL = 0.d0/SQRT(FOURPI)
    byL = 1.d0 ! 3.d0/SQRT(FOURPI)
    bzL = 0.d0/SQRT(FOURPI)
    
    ! Riemann problem right state
    rhoR = 1.d0
    vxR = 1.d0
    vyR = 0.d0
    vzR = 0.d0
    bxR = bxL
    byR = 1.d0 ! 0.d0/SQRT(FOURPI)
    bzR = 0.d0/SQRT(FOURPI)
    
    DO i = 1-nb, nx+nb
         
            
        IF(i .LT. nx/2) THEN
        
            q(i,1) = rhoL
            q(i,2) = rhoL*vxL        
            q(i,3) = rhoL*vyL
            q(i,4) = rhoL*vzL    
            q(i,5) = bxL    
            q(i,6) = byL    
            q(i,7) = bzL     
            
        ELSE
            
            q(i,1) = rhoR
            q(i,2) = rhoR*vxR        
            q(i,3) = rhoR*vyR
            q(i,4) = rhoR*vzR    
            q(i,5) = bxR    
            q(i,6) = byR    
            q(i,7) = bzR     
            
        END IF 
                
    END DO
    
    CALL fileOutput()

    
END SUBROUTINE init_state



SUBROUTINE fileOutput()

    INTEGER :: i
    REAL*8 :: dx, vx, vy, vz
    
    dx = 1.d0/nx

    DO i =1,nx
        vx = q(i,2)/q(i,1)
        vy = q(i,3)/q(i,1)
        vz = q(i,4)/q(i,1)
        WRITE(10,*) i*dx,q(i,1),vx,vy,vz
        WRITE(11,*) I*dx,q(i,6),q(i,7)        
    END DO


END SUBROUTINE fileOutput


SUBROUTINE boundary_conditions()

     INTEGER :: i
     
    ! Open Boundaries
    DO i=1,nb
        q(1-i,:) = q(1,:)       
        q(nx+i,:) = q(nx,:)
    END DO


END SUBROUTINE boundary_conditions


END PROGRAM isothermal_tvd_driver