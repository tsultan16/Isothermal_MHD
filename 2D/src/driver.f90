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
OPEN(UNIT=14, FILE='Output/fluid_diagcut.txt')
OPEN(UNIT=15, FILE='Output/Bfield_diagcut.txt')


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
    REAL*8 :: dx, x, y, isq2
    REAL*8 :: rhoL, rhoR, vxL, vxR, vyL, vyR, vzL, vzR, &
              bxL, bxR, byL, byR, bzL, bzR
    REAL*8 :: vL_par, vL_perp, BL_par, BL_perp, vR_par, vR_perp, BR_par, BR_perp
    
    INTEGER, PARAMETER :: dir = 3 ! 1 : x, 2 : y, 3 : diagonal  
    
    
    dx = 1.d0 / nx
    
    IF(dir .EQ. 1) THEN
        
    ! Riemann problem left state
    rhoL = 1.d0
    vxL = -1.0d0
    vyL = 0.d0
    vzL = 0.d0
    bxL = 0.d0/SQRT(FOURPI)
    byL = 1.d0 ! 3.6d0/SQRT(FOURPI)
    bzL = 0.0d0/SQRT(FOURPI)
    
    ! Riemann problem right state
    rhoR = 1.d0
    vxR = 1.d0
    vyR = 0.d0
    vzR = 0.d0
    bxR = bxL
    byR = 1.d0 ! 1.d0/SQRT(FOURPI)
    bzR = 0.d0/SQRT(FOURPI)
    
    
    DO j = 1-nb, ny+nb
    DO i = 1-nb, nx+nb
                   
        IF(i .LT. nx/2) THEN
        
            q_2d(i,j,1) = rhoL
            q_2d(i,j,2) = rhoL*vxL        
            q_2d(i,j,3) = rhoL*vyL
            q_2d(i,j,4) = rhoL*vzL    
            q_2d(i,j,7) = bzL     
  
!            q_2d(i,j,5) = bxL
!            q_2d(i,j,6) = byL
  
  
            bface_2d(i,j,1) = bxL
            bface_2d(i,j,2) = byL
            
        ELSE
            
            q_2d(i,j,1) = rhoR
            q_2d(i,j,2) = rhoR*vxR        
            q_2d(i,j,3) = rhoR*vyR
            q_2d(i,j,4) = rhoR*vzR    
            q_2d(i,j,7) = bzR     
           
            bface_2d(i,j,1) = bxR
            bface_2d(i,j,2) = byR

!            q_2d(i,j,5) = bxR
!            q_2d(i,j,6) = byR
           
        END IF 
                
    END DO
    END DO

    
    END IF    
    
    
    IF(dir .EQ. 2) THEN
    
    ! Riemann problem left state
    rhoL = 1.d0
    vyL = -1.0d0
    vxL = 0.d0
    vzL = 0.d0
    byL = 0.d0/SQRT(FOURPI)
    bxL = 1.d0 ! 3.6d0/SQRT(FOURPI)
    bzL = 0.0d0/SQRT(FOURPI)
    
    ! Riemann problem right state
    rhoR = 1.d0
    vyR = 1.d0
    vxR = 0.d0
    vzR = 0.d0
    byR = byL
    bxR = 1.d0 ! 1.d0/SQRT(FOURPI)
    bzR = 0.d0/SQRT(FOURPI)
    
    
    DO j = 1-nb, ny+nb
    DO i = 1-nb, nx+nb
                   
        IF(j .LT. ny/2) THEN
        
            q_2d(i,j,1) = rhoL
            q_2d(i,j,2) = rhoL*vxL        
            q_2d(i,j,3) = rhoL*vyL
            q_2d(i,j,4) = rhoL*vzL    
            q_2d(i,j,7) = bzL     
  
!            q_2d(i,j,5) = bxL
!            q_2d(i,j,6) = byL
  
  
            bface_2d(i,j,1) = bxL
            bface_2d(i,j,2) = byL
            
        ELSE
            
            q_2d(i,j,1) = rhoR
            q_2d(i,j,2) = rhoR*vxR        
            q_2d(i,j,3) = rhoR*vyR
            q_2d(i,j,4) = rhoR*vzR    
            q_2d(i,j,7) = bzR     

!            q_2d(i,j,5) = bxR
!            q_2d(i,j,6) = byR
            
            bface_2d(i,j,1) = bxR
            bface_2d(i,j,2) = byR
          
           
        END IF 
                
    END DO
    END DO    
    
    END IF  
    


    IF(dir .EQ. 3) THEN
    
    ! Riemann problem left state
    rhoL = 1.08d0
    vL_perp = 1.2d0
    vL_par = 0.01d0
    vzL = 0.5d0
    bL_perp = 2.d0/SQRT(FOURPI)
    bL_par = 3.6d0/SQRT(FOURPI)
    bzL = 2.0d0/SQRT(FOURPI)
    
    ! Riemann problem right state
    rhoR = 1.0d0
    vR_perp = 0.d0
    vR_par = 0.d0
    vzR = 0.d0
    bR_perp = bL_perp
    bR_par = 4.d0/SQRT(FOURPI)
    bzR = 2.d0/SQRT(FOURPI)
    
    isq2 = 1.d0 /SQRT(2.d0)
    
    DO j = 1-nb, ny+nb
    DO i = 1-nb, nx+nb
               
        x = i * dx
        y = j * dx
               
        IF(y .LT. 1.d0 - x) THEN
        
            q_2d(i,j,1) = rhoL
            q_2d(i,j,2) = rhoL * isq2 *(vL_perp - vL_par)        
            q_2d(i,j,3) = rhoL * isq2 *(vL_perp + vL_par)
            q_2d(i,j,4) = rhoL * vzL    
            q_2d(i,j,7) = bzL     

            bface_2d(i,j,1) = isq2 *(bL_perp - bL_par)  
            bface_2d(i,j,2) = isq2 *(bL_perp + bL_par)  
            
        ELSE
            
            q_2d(i,j,1) = rhoR
            q_2d(i,j,2) = rhoR * isq2 *(vR_perp - vR_par)        
            q_2d(i,j,3) = rhoR * isq2 *(vR_perp + vR_par)
            q_2d(i,j,4) = rhoR * vzR    
            q_2d(i,j,7) = bzR     

!            q_2d(i,j,5) = bxR
!            q_2d(i,j,6) = byR
            
            bface_2d(i,j,1) = isq2 *(bR_perp - bR_par)  
            bface_2d(i,j,2) = isq2 *(bR_perp + bR_par)  
          
           
        END IF 
                
    END DO
    END DO    
    
    END IF  

    
    DO j = 1-nb+1, ny+nb    
    DO i = 1-nb+1, nx+nb

        q_2d(i,j,5) = 0.5d0 * ( bface_2d(i-1,j,1) + bface_2d(i,j,1) )     
        q_2d(i,j,6) = 0.5d0 * ( bface_2d(i,j-1,2) + bface_2d(i,j,2) )       
        
    END DO
    END DO
    
    
    
    CALL boundary_conditions()
    
    CALL fileOutput()

    
END SUBROUTINE init_state


SUBROUTINE fileOutput()

    INTEGER :: i, j
    REAL*8 :: dx, vx, vy, vz, vpar, vperp, bpar, bperp, isqr2
    
    dx = 1.d0/nx
    j = ny/2
    
    isqr2 = 1.d0 / SQRT(2.d0)
    
    
    DO i =1,nx
        vx = q_2d(i,j,2)/q_2d(i,j,1)
        vy = q_2d(i,j,3)/q_2d(i,j,1)
        vz = q_2d(i,j,4)/q_2d(i,j,1)
        WRITE(10,*) i*dx,q_2d(i,j,1),vx,vy,vz
        WRITE(11,*) i*dx,q_2d(i,j,5),q_2d(i,j,6),q_2d(i,j,7)        
    END DO

    i = nx/2
    
    DO j =1,ny
        vx = q_2d(i,j,2)/q_2d(i,j,1)
        vy = q_2d(i,j,3)/q_2d(i,j,1)
        vz = q_2d(i,j,4)/q_2d(i,j,1)
        WRITE(12,*) j*dx,q_2d(i,j,1),vx,vy,vz
        WRITE(13,*) j*dx,q_2d(i,j,5),q_2d(i,j,6),q_2d(i,j,7)        
    END DO


    DO i =1,nx
        vx = q_2d(i,i,2)/q_2d(i,i,1)
        vy = q_2d(i,i,3)/q_2d(i,i,1)
        vz = q_2d(i,i,4)/q_2d(i,i,1)
        
        vpar = isqr2*(vy-vx)
        vperp = isqr2*(vy+vx)
        bpar = isqr2*(q_2d(i,i,6) - q_2d(i,i,5))
        bperp = isqr2*(q_2d(i,i,6) + q_2d(i,i,5))
        
        WRITE(14,*) i*SQRT(2.d0)*dx,q_2d(i,i,1),vperp,vpar,vz
        WRITE(15,*) i*SQRT(2.d0)*dx,bperp,bpar,q_2d(i,i,7)        
    END DO



END SUBROUTINE fileOutput


SUBROUTINE boundary_conditions()

     INTEGER :: i, j, k
     
    ! Open Boundaries
    
    DO j=1-nb,ny+nb
            
        q_2d(1-nb+2,j,:) = q_2d(1-nb+3,j,:)     
        q_2d(1-nb+1,j,:) = q_2d(1-nb+3,j,:)     
        q_2d(1-nb,j,:)   = q_2d(1-nb+3,j,:)
        
        bface_2d(1-nb+2,j,:) = bface_2d(1-nb+3,j,:)     
        bface_2d(1-nb+1,j,:) = bface_2d(1-nb+3,j,:)     
        bface_2d(1-nb,j,:)   = bface_2d(1-nb+3,j,:)        
        
        q_2d(nx+nb-2,j,:) = q_2d(nx+nb-3,j,:)
        q_2d(nx+nb-1,j,:) = q_2d(nx+nb-3,j,:)
        q_2d(nx+nb,j,:)   = q_2d(nx+nb-3,j,:)
        
        bface_2d(nx+nb-2,j,:) = bface_2d(nx+nb-3,j,:)
        bface_2d(nx+nb-1,j,:) = bface_2d(nx+nb-3,j,:)
        bface_2d(nx+nb,j,:)   = bface_2d(nx+nb-3,j,:)

    END DO

    DO i=1-nb,nx+nb

        q_2d(i,1-nb+2,:) = q_2d(i,1-nb+3,:)       
        q_2d(i,1-nb+1,:) = q_2d(i,1-nb+3,:)       
        q_2d(i,1-nb,:)   = q_2d(i,1-nb+3,:)
        
        bface_2d(i,1-nb+2,:) = bface_2d(i,1-nb+3,:)       
        bface_2d(i,1-nb+1,:) = bface_2d(i,1-nb+3,:)       
        bface_2d(i,1-nb,:)   = bface_2d(i,1-nb+3,:)     
      
        q_2d(i,ny+nb-2,:) = q_2d(i,ny+nb-3,:)
        q_2d(i,ny+nb-1,:) = q_2d(i,ny+nb-3,:)
        q_2d(i,ny+nb,:)   = q_2d(i,ny+nb-3,:)

        bface_2d(i,ny+nb-2,:) = bface_2d(i,ny+nb-3,:)
        bface_2d(i,ny+nb-1,:) = bface_2d(i,ny+nb-3,:)
        bface_2d(i,ny+nb,:)   = bface_2d(i,ny+nb-3,:)

    END DO

  


    


END SUBROUTINE boundary_conditions


END PROGRAM isothermal_tvd_driver