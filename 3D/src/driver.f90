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
OPEN(UNIT=14, FILE='Output/fluid_zcut.txt')
OPEN(UNIT=15, FILE='Output/Bfield_zcut.txt')
OPEN(UNIT=16, FILE='Output/fluid_xydiagcut.txt')
OPEN(UNIT=17, FILE='Output/Bfield_xydiagcut.txt')
OPEN(UNIT=18, FILE='Output/fluid_yzdiagcut.txt')
OPEN(UNIT=19, FILE='Output/Bfield_yzdiagcut.txt')

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
    CALL solve_3d()
    
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
CLOSE(UNIT=12)
CLOSE(UNIT=13)
CLOSE(UNIT=14)
CLOSE(UNIT=15)
CLOSE(UNIT=16)
CLOSE(UNIT=17)
CLOSE(UNIT=18)
CLOSE(UNIT=19)

PRINT*,'Done!'

CONTAINS


SUBROUTINE init_state()

    INTEGER :: i, j, k
    REAL*8 :: dx, x, y, z, isq2, tmp
    REAL*8 :: rhoL, rhoR, vxL, vxR, vyL, vyR, vzL, vzR, &
              bxL, bxR, byL, byR, bzL, bzR
    REAL*8 :: vL_par, vL_perp, BL_par, BL_perp, vR_par, vR_perp, BR_par, BR_perp
    
    INTEGER, PARAMETER :: dir = 5 ! 1 : x , 2 : y , 3 : z , 4: x=y , 5: y=z  
    INTEGER, PARAMETER :: test = 1
    
    dx = 1.d0 / MAX(nx,ny,nz)
    
    
    
    IF(test .EQ. 1) THEN
        
        ! Riemann problem left state
        rhoL = 1.08d0
        vxL  = 1.2d0
        vyL  = 0.01d0
        vzL  = 0.5d0
        bxL  = 2.d0/SQRT(FOURPI)
        byL  = 3.6d0/SQRT(FOURPI)
        bzL  = 2.d0/SQRT(FOURPI)
    
        ! Riemann problem right state
        rhoR = 1.d0
        vxR  = 0.d0
        vyR  = 0.d0
        vzR  = 0.d0
        bxR  = bxL
        byR  = 4.d0/SQRT(FOURPI)
        bzR  = 2.d0/SQRT(FOURPI)
    
    
    ELSE IF(test .EQ. 2) THEN
    
        ! Riemann problem left state
        rhoL = 1.d0
        vxL  = -1.d0
        vyL  = 0.d0
        vzL  = 0.d0
        bxL  = 0.d0
        byL  = 1.d0
        bzL  = 0.d0
    
        ! Riemann problem right state
        rhoR = 1.d0
        vxR  = 1.d0
        vyR  = 0.d0
        vzR  = 0.d0
        bxR  = bxL
        byR  = 1.d0
        bzR  = 0.d0
    
    END IF
    
    
    IF(dir .EQ. 1) THEN
        
    DO k = 1-nb, nz+nb
    DO j = 1-nb, ny+nb
    DO i = 1-nb, nx+nb
                   
        IF(i .LT. nx/2) THEN
        
            q_3d(i,j,k,1) = rhoL
            q_3d(i,j,k,2) = rhoL*vxL        
            q_3d(i,j,k,3) = rhoL*vyL
            q_3d(i,j,k,4) = rhoL*vzL     
            bface_3d(i,j,k,1) = bxL
            bface_3d(i,j,k,2) = byL
            bface_3d(i,j,k,3) = bzL
            
        ELSE
            
            q_3d(i,j,k,1) = rhoR
            q_3d(i,j,k,2) = rhoR*vxR        
            q_3d(i,j,k,3) = rhoR*vyR
            q_3d(i,j,k,4) = rhoR*vzR    
            bface_3d(i,j,k,1) = bxR
            bface_3d(i,j,k,2) = byR
            bface_3d(i,j,k,3) = bzR


        END IF 
                
    END DO
    END DO
    END DO

    
    END IF    
    
    
    IF(dir .EQ. 2) THEN
    
    tmp = vxL
    vxL = vyL
    vyL = tmp  
    
    tmp = BxL
    BxL = ByL
    ByL = tmp

    tmp = vxR
    vxR = vyR
    vyR = tmp  
    
    tmp = BxR
    BxR = ByR
    ByR = tmp
    
    
    DO k = 1-nb, nz+nb
    DO j = 1-nb, ny+nb
    DO i = 1-nb, nx+nb
                   
        IF(j .LT. ny/2) THEN
        
            q_3d(i,j,k,1) = rhoL
            q_3d(i,j,k,2) = rhoL*vxL        
            q_3d(i,j,k,3) = rhoL*vyL
            q_3d(i,j,k,4) = rhoL*vzL     
            bface_3d(i,j,k,1) = bxL
            bface_3d(i,j,k,2) = byL
            bface_3d(i,j,k,3) = bzL
            
            
        ELSE
            
            q_3d(i,j,k,1) = rhoR
            q_3d(i,j,k,2) = rhoR*vxR        
            q_3d(i,j,k,3) = rhoR*vyR
            q_3d(i,j,k,4) = rhoR*vzR    
            bface_3d(i,j,k,1) = bxR
            bface_3d(i,j,k,2) = byR
            bface_3d(i,j,k,3) = bzR
          
           
        END IF 
                
    END DO
    END DO    
    END DO    
    
    END IF  
    
    
     IF(dir .EQ. 3) THEN
    
    tmp = vxL
    vxL = vzL
    vzL = tmp  
    
    tmp = BxL
    BxL = BzL
    BzL = tmp

    tmp = vxR
    vxR = vzR
    vzR = tmp  
    
    tmp = BxR
    BxR = BzR
    BzR = tmp
    
    
    DO k = 1-nb, nz+nb
    DO j = 1-nb, ny+nb
    DO i = 1-nb, nx+nb
                   
        IF(k .LT. nz/2) THEN
        
            q_3d(i,j,k,1) = rhoL
            q_3d(i,j,k,2) = rhoL*vxL        
            q_3d(i,j,k,3) = rhoL*vyL
            q_3d(i,j,k,4) = rhoL*vzL     
            bface_3d(i,j,k,1) = bxL
            bface_3d(i,j,k,2) = byL
            bface_3d(i,j,k,3) = bzL
            
            
        ELSE
            
            q_3d(i,j,k,1) = rhoR
            q_3d(i,j,k,2) = rhoR*vxR        
            q_3d(i,j,k,3) = rhoR*vyR
            q_3d(i,j,k,4) = rhoR*vzR    
            bface_3d(i,j,k,1) = bxR
            bface_3d(i,j,k,2) = byR
            bface_3d(i,j,k,3) = bzR
          
           
        END IF 
                
    END DO
    END DO    
    END DO    
    
    END IF  
    
    
    IF(dir .EQ. 4) THEN
   
        vL_perp = vxL
        vL_par  = vyL
        bL_perp = bxL
        bL_par  = byL
   
        vR_perp = vxR
        vR_par  = vyR
        bR_perp = bL_perp
        bR_par  = byR
    
    
        isq2 = 1.d0 /SQRT(2.d0)
    
        DO k = 1-nb, nz+nb
        DO j = 1-nb, ny+nb
       
        y = j * dx
        
        DO i = 1-nb, nx+nb
               
        x = i * dx
               
        IF(y .LT. 1.d0 - x) THEN
        
            q_3d(i,j,k,1) = rhoL
            q_3d(i,j,k,2) = rhoL * isq2 *(vL_perp - vL_par)        
            q_3d(i,j,k,3) = rhoL * isq2 *(vL_perp + vL_par)
            q_3d(i,j,k,4) = rhoL * vzL    
            
            bface_3d(i,j,k,1) = isq2 *(bL_perp - bL_par)  
            bface_3d(i,j,k,2) = isq2 *(bL_perp + bL_par)  
            bface_3d(i,j,k,3) = bzL  
            
        ELSE
            
            q_3d(i,j,k,1) = rhoR
            q_3d(i,j,k,2) = rhoR * isq2 *(vR_perp - vR_par)        
            q_3d(i,j,k,3) = rhoR * isq2 *(vR_perp + vR_par)
            q_3d(i,j,k,4) = rhoR * vzR    


            bface_3d(i,j,k,1) = isq2 *(bR_perp - bR_par)  
            bface_3d(i,j,k,2) = isq2 *(bR_perp + bR_par)  
            bface_3d(i,j,k,3) = bzR
          
           
        END IF 
                
    END DO
    END DO    
    END DO    
    
    END IF  


    IF(dir .EQ. 5) THEN
       
        vL_perp = vxL
        vL_par  = vyL
        bL_perp = bxL
        bL_par  = byL
   
        vR_perp = vxR
        vR_par  = vyR
        bR_perp = bL_perp
        bR_par  = byR
    
    
        isq2 = 1.d0 /SQRT(2.d0)
    
        DO k = 1-nb, nz+nb
        
        z = k * dx
        
        DO j = 1-nb, ny+nb
       
        y = j * dx
        
        DO i = 1-nb, nx+nb
               
        IF(z .LT. 1.d0 - y) THEN
        
            q_3d(i,j,k,1) = rhoL
            q_3d(i,j,k,3) = rhoL * isq2 *(vL_perp - vL_par)        
            q_3d(i,j,k,4) = rhoL * isq2 *(vL_perp + vL_par)
            q_3d(i,j,k,2) = rhoL * vzL    
            
            bface_3d(i,j,k,2) = isq2 *(bL_perp - bL_par)  
            bface_3d(i,j,k,3) = isq2 *(bL_perp + bL_par)  
            bface_3d(i,j,k,1) = bzL  
            
        ELSE
            
            q_3d(i,j,k,1) = rhoR
            q_3d(i,j,k,3) = rhoR * isq2 *(vR_perp - vR_par)        
            q_3d(i,j,k,4) = rhoR * isq2 *(vR_perp + vR_par)
            q_3d(i,j,k,2) = rhoR * vzR    


            bface_3d(i,j,k,2) = isq2 *(bR_perp - bR_par)  
            bface_3d(i,j,k,3) = isq2 *(bR_perp + bR_par)  
            bface_3d(i,j,k,1) = bzR
          
           
        END IF 
                
        END DO
        END DO    
        END DO    
    
    END IF
    
    
    DO k = 1-nb+1, nz+nb    
    DO j = 1-nb+1, ny+nb    
    DO i = 1-nb+1, nx+nb

        q_3d(i,j,k,5) = 0.5d0 * ( bface_3d(i-1,j,k,1) + bface_3d(i,j,k,1) )     
        q_3d(i,j,k,6) = 0.5d0 * ( bface_3d(i,j-1,k,2) + bface_3d(i,j,k,2) )       
        q_3d(i,j,k,7) = 0.5d0 * ( bface_3d(i,j,k-1,3) + bface_3d(i,j,k,3) )       
        
    END DO
    END DO
    END DO
    
    
    
    CALL boundary_conditions()
    
    CALL fileOutput()

    
END SUBROUTINE init_state


SUBROUTINE fileOutput()

    INTEGER :: i, j, k
    REAL*8 :: dx, vx, vy, vz, vpar, vperp, bpar, bperp, isqr2
    
    dx = 1.d0 / MAX(nx,ny,nz)
    isqr2 = 1.d0 / SQRT(2.d0)
    
    j = ny/2
    k = nz/2
    
    DO i =1,nx
        vx = q_3d(i,j,k,2)/q_3d(i,j,k,1)
        vy = q_3d(i,j,k,3)/q_3d(i,j,k,1)
        vz = q_3d(i,j,k,4)/q_3d(i,j,k,1)
        WRITE(10,*) i*dx,q_3d(i,j,k,1),vx,vy,vz
        WRITE(11,*) i*dx,q_3d(i,j,k,5),q_3d(i,j,k,6),q_3d(i,j,k,7)        
    END DO

    i = nx/2
    k = nz/2
    
    DO j =1,ny
        vx = q_3d(i,j,k,2)/q_3d(i,j,k,1)
        vy = q_3d(i,j,k,3)/q_3d(i,j,k,1)
        vz = q_3d(i,j,k,4)/q_3d(i,j,k,1)
        WRITE(12,*) j*dx,q_3d(i,j,k,1),vx,vy,vz
        WRITE(13,*) j*dx,q_3d(i,j,k,5),q_3d(i,j,k,6),q_3d(i,j,k,7)        
    END DO

    i = nx/2
    j = ny/2
    
    DO k =1,nz
        vx = q_3d(i,j,k,2)/q_3d(i,j,k,1)
        vy = q_3d(i,j,k,3)/q_3d(i,j,k,1)
        vz = q_3d(i,j,k,4)/q_3d(i,j,k,1)
        WRITE(14,*) k*dx,q_3d(i,j,k,1),vx,vy,vz
        WRITE(15,*) k*dx,q_3d(i,j,k,5),q_3d(i,j,k,6),q_3d(i,j,k,7)        
    END DO

    k = nz/2
    
    DO j =1,ny
    DO i =1,nx
        
        IF(i .EQ. j) THEN 

            vx = q_3d(i,j,k,2)/q_3d(i,j,k,1)
            vy = q_3d(i,j,k,3)/q_3d(i,j,k,1)
            vz = q_3d(i,j,k,4)/q_3d(i,j,k,1)

            vperp = isqr2 * (vy + vx) 
            vpar  = isqr2 * (vy - vx)
            bperp = isqr2 * (q_3d(i,j,k,6) + q_3d(i,j,k,5))
            bpar  = isqr2 * (q_3d(i,j,k,6) - q_3d(i,j,k,5))            
        
            WRITE(16,*) i*dx,q_3d(i,j,k,1),vperp,vpar,vz
            WRITE(17,*) i*dx,bperp,bpar,q_3d(i,j,k,7)
            
        END IF      
        
    END DO
    END DO



    i = nx/2
    
    DO k =1,nz
    DO j =1,ny
        
        IF(j .EQ. k) THEN 

            vx = q_3d(i,j,k,2)/q_3d(i,j,k,1)
            vy = q_3d(i,j,k,3)/q_3d(i,j,k,1)
            vz = q_3d(i,j,k,4)/q_3d(i,j,k,1)

            vperp = isqr2 * (vz + vy) 
            vpar  = isqr2 * (vz - vy)
            bperp = isqr2 * (q_3d(i,j,k,7) + q_3d(i,j,k,6))
            bpar  = isqr2 * (q_3d(i,j,k,7) - q_3d(i,j,k,6))            
        
            WRITE(18,*) j*dx,q_3d(i,j,k,1),vperp,vpar,vx
            WRITE(19,*) j*dx,bperp,bpar,q_3d(i,j,k,5)
            
        END IF      
        
    END DO
    END DO
    

END SUBROUTINE fileOutput


! apply boundary conditions to state vector and cell corner electric fields
SUBROUTINE boundary_conditions()

     INTEGER :: i, j, k
     
     ! Open Boundaries
    
        !************        
        ! x boundary
        !************        
        
        DO i = 0, nb-1
            q_3d(1-nb+i,:,:,:) = q_3d(1,:,:,:)     
            q_3d(nx+1+i,:,:,:) = q_3d(nx,:,:,:)      

            bface_3d(1-nb+i,:,:,:) = bface_3d(1,:,:,:)     
            bface_3d(nx+1+i,:,:,:) = bface_3d(nx,:,:,:)      

            !emf_corner(1-nb+i,:,:,:) = emf_corner(1,:,:,:)     
            !emf_corner(nx+1+i,:,:,:) = emf_corner(nx,:,:,:)      

        END DO
        
  
        !************
        ! y boundary
        !************

        DO j = 0, nb-1

            q_3d(:,1-nb+j,:,:) = q_3d(:,1,:,:)     
            q_3d(:,ny+1+j,:,:) = q_3d(:,ny,:,:)    
            
            bface_3d(:,1-nb+j,:,:) = bface_3d(:,1,:,:)     
            bface_3d(:,ny+1+j,:,:) = bface_3d(:,ny,:,:)    
            
            !emf_corner(:,1-nb+j,:,:) = emf_corner(:,1,:,:)     
            !emf_corner(:,ny+1+j,:,:) = emf_corner(:,ny,:,:)    
            
        END DO      
        

        !************ 
        ! z boundary
        !************

        DO k = 0, nb-1

            q_3d(:,:,1-nb+k,:) = q_3d(:,:,1,:)     
            q_3d(:,:,nz+1+k,:) = q_3d(:,:,nz,:)      

            bface_3d(:,:,1-nb+k,:) = bface_3d(:,:,1,:)     
            bface_3d(:,:,nz+1+k,:) = bface_3d(:,:,nz,:)      
            
            !emf_corner(:,:,1-nb+k,:) = emf_corner(:,:,1,:)     
            !emf_corner(:,:,nz+1+k,:) = emf_corner(:,:,nz,:)  
            
            
        END DO
        
     
        

END SUBROUTINE boundary_conditions


END PROGRAM isothermal_tvd_driver