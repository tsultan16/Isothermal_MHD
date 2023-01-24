MODULE init_domain_mod

USE constants_mod
USE grid_data_mod

IMPLICIT NONE

CONTAINS




SUBROUTINE init_state()

    INTEGER :: i, j, k
    REAL*8 :: dx, x, y, z, isq2, tmp
    REAL*8 :: rhoL, rhoR, vxL, vxR, vyL, vyR, vzL, vzR, &
              bxL, bxR, byL, byR, bzL, bzR
    REAL*8 :: vL_par, vL_perp, BL_par, BL_perp, vR_par, vR_perp, BR_par, BR_perp
    
    INTEGER, PARAMETER :: dir = 1 ! 1 : x , 2 : y , 3 : z , 4: x=y , 5: y=z  
    INTEGER, PARAMETER :: test = 0
    
    dx = 1.d0 / MAX(nx,ny,nz)
    
    
    IF(test .EQ. 0) THEN
    
        ! Riemann problem left state
        rhoL = 1.d0
        vxL  = 1.d0
        vyL  = 0.d0
        vzL  = 0.d0
        bxL  = 0.d0
        byL  = 0.d0
        bzL  = 0.d0
    
        ! Riemann problem right state
        rhoR = 1.d0
        vxR  = 1.d0
        vyR  = 0.d0
        vzR  = 0.d0
        bxR  = 0.d0
        byR  = 0.d0
        bzR  = 0.d0    
    
           
    END IF
    
    
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
        
    
         !IF(my_coord(1) .EQ. 0) THEN
       
        !IF(i .LT. nx/2) THEN
        
            q_3d(i,j,k,1) = rhoL
            q_3d(i,j,k,2) = rhoL*vxL        
            q_3d(i,j,k,3) = rhoL*vyL
            q_3d(i,j,k,4) = rhoL*vzL     
            q_3d(i,j,k,5) = bxL    
            q_3d(i,j,k,6) = byL     
            q_3d(i,j,k,7) = bzL
            
            bface_3d(i,j,k,1) = bxL
            bface_3d(i,j,k,2) = byL
            bface_3d(i,j,k,3) = bzL
          
          !END IF
          
          
          !IF(my_coord(1) .EQ. 1)THEN
            
          !  IF(i .LT. nx/2) THEN

            !q_3d(i,j,k,1) = rhoL
            !q_3d(i,j,k,2) = rhoL*vxL        
            !q_3d(i,j,k,3) = rhoL*vyL
            !q_3d(i,j,k,4) = rhoL*vzL     
            !q_3d(i,j,k,5) = bxL    
            !q_3d(i,j,k,6) = byL     
            !q_3d(i,j,k,7) = bzL
            
            !bface_3d(i,j,k,1) = bxL
            !bface_3d(i,j,k,2) = byL
            !bface_3d(i,j,k,3) = bzL
            
            !ELSE            
            
            !q_3d(i,j,k,1) = rhoR
            !q_3d(i,j,k,2) = rhoR*vxR        
            !q_3d(i,j,k,3) = rhoR*vyR
            !q_3d(i,j,k,4) = rhoR*vzR    
            !q_3d(i,j,k,5) = bxR    
            !q_3d(i,j,k,6) = byR     
            !q_3d(i,j,k,7) = bzR
            
            !bface_3d(i,j,k,1) = bxR
            !bface_3d(i,j,k,2) = byR
            !bface_3d(i,j,k,3) = bzR
            
            !END IF

          !END IF
  
          !IF(my_coord(1) .EQ. 2)THEN
            
            !q_3d(i,j,k,1) = rhoR
            !q_3d(i,j,k,2) = rhoR*vxR        
            !q_3d(i,j,k,3) = rhoR*vyR
            !q_3d(i,j,k,4) = rhoR*vzR    
            !q_3d(i,j,k,5) = bxR    
            !q_3d(i,j,k,6) = byR     
            !q_3d(i,j,k,7) = bzR
            
            !bface_3d(i,j,k,1) = bxR
            !bface_3d(i,j,k,2) = byR
            !bface_3d(i,j,k,3) = bzR

          !END IF      
        
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
   
    
    
    
END SUBROUTINE init_state


END MODULE init_domain_mod