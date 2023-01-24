! 1D Isothermal MHD 2nd order TVD Solver
! References: Kim 'et al, 1998  (K98)


MODULE Iso_TVD_Solver_mod

USE constants_mod
USE grid_data_mod

IMPLICIT NONE

! Solver Parameters
INTEGER, PARAMETER :: nwaves = 6
REAL*8, PARAMETER  :: ZERO = 1.D-20 
REAL*8, PARAMETER :: limiter_beta = 1.D0  ! Sweby limiter: 1 <= beta <= 2
REAL*8, PARAMETER :: COUR = 0.8D0 ! Courant number

! Solver Local variables 
REAL*8, ALLOCATABLE :: Fs(:,:),    &        ! time-averaged flux vector
                       qavg(:,:),  &        ! Roe-averaged state    
                       Rk(:,:,:),  &        ! Roe-averaged right eigenvectors
                       ck(:,:),    &        ! characteristic variables (i.e. projection of state vector onto the left eigenvectors)   
                       eigenvalues(:,:),  & ! Roe-averaged eigenvalues
                       wave_speeds(:,:)     ! Roe-averaged wave speeds
                       
REAL*8, PARAMETER :: eps(6) = (/ 0.3, 0.0, 0.3, 0.3, 0.0 ,0.3/) ! wave dissipation constants

REAL*8 :: dt, dx, dtdx


CONTAINS


SUBROUTINE init_solver()
    
    OPEN(UNIT = 12, FILE='Output/flux.txt')
    CALL allocate_local_variables()


END SUBROUTINE init_solver


SUBROUTINE destroy_solver()

   CALL deallocate_local_variables()
   CLOSE(UNIT=12)
   
END SUBROUTINE destroy_solver

! allocate local variables/work arrays
SUBROUTINE allocate_local_variables()

    ALLOCATE(Fs(0:nx,nwaves))
    ALLOCATE(qavg(1-nb:nx+nb, 7))
    ALLOCATE(Rk(1-nb:nx+nb, nwaves, nwaves))
    ALLOCATE(ck(1-nb:nx+nb, nwaves))
    ALLOCATE(eigenvalues(1-nb:nx+nb, nwaves))
    ALLOCATE(wave_speeds(1-nb:nx+nb, 3))

END SUBROUTINE allocate_local_variables


SUBROUTINE deallocate_local_variables()

    DEALLOCATE(Fs, qavg, Rk, ck, eigenvalues, wave_speeds)
    
END SUBROUTINE deallocate_local_variables

! Top-level Isothermal TVD MHD solve routine
SUBROUTINE solve()

    REAL*8 :: cmax

    dx = 1.d0/nx
    dt = 0.d0
    
    ! compute Roe-averaged state at cell-interface
    PRINT*,'Computing Roe averaged variables.'
    CALL compute_avg_state()
    
    ! compute Roe-averaged eigenvalues
    PRINT*,'Computing eigenvalues.'
    CALL compute_eigenvalues()
    
    ! compute time-step size
    PRINT*,'Computing time step size.'
    CALL get_max_wavespeed(cmax)
    
    dt = COUR * dx / cmax 
    dtdx = dt/dx
    PRINT*,'dt =',dt

    ! compute Roe-averaged right eigenvectors and characteristic variables
    PRINT*,'Computing eigenvectors and characteristics.'
    CALL compute_eigenvectors()

    ! compute time-averaged cell-interface (2nd order) TVD fluxes
    PRINT*,'Computing TVD fluxes'
    CALL compute_TVD_fluxes()
 
    ! update state vector
    PRINT*,'Updating state vector.'
    CALL update_state_vector()


END SUBROUTINE solve



! Compute cell-interface roe averaged state
SUBROUTINE compute_avg_state()

    REAL*8 :: rho2, vx2, vy2, vz2, Bx2, By2, Bz2, isrho
    INTEGER :: i
    
    
    DO i = 1-nb, nx+nb-1

        ! compute Roe-averaged primitive variables
        rho2 = 0.5d0 * ( q(i,1) + q(i+1,1) )
        isrho = 1.D0 / SQRT(rho2)
        vx2  = 0.5d0 * ( q(i,2)/q(i,1) + q(i+1,2)/q(i+1,1) )
        vy2  = 0.5d0 * ( q(i,3)/q(i,1) + q(i+1,3)/q(i+1,1) )
        vz2  = 0.5d0 * ( q(i,4)/q(i,1) + q(i+1,4)/q(i+1,1) )
        Bx2  = 0.5d0 * ( q(i,5)+ q(i+1,5) )
        By2  = 0.5d0 * ( q(i,6)+ q(i+1,6) )
        Bz2  = 0.5d0 * ( q(i,7)+ q(i+1,7) )
        
        qavg(i,1) = rho2
        qavg(i,2) = vx2
        qavg(i,3) = vy2
        qavg(i,4) = vz2
        qavg(i,5) = Bx2 * isrho
        qavg(i,6) = By2 * isrho
        qavg(i,7) = Bz2 * isrho           

    END DO 



END SUBROUTINE compute_avg_state


! Compute eigenvalues
SUBROUTINE compute_eigenvalues()

    REAL*8 :: rho2, vx2, vy2, vz2, bx2, by2, bz2, asqr
    REAL*8 :: cf, ca, cs, bsqr, tmp1, tmp2, isrho
    INTEGER :: i
        
    DO i = 1-nb, nx+nb-1

        rho2 = qavg(i,1)
        vx2  = qavg(i,2)
        asqr = sound_speed**2
        bsqr = qavg(i,5)**2 + qavg(i,6)**2 + qavg(i,7)**2  
        tmp1 = asqr + bsqr        
        tmp2 = SQRT(tmp1**2 - 4.d0 * asqr * qavg(i,5)**2 )

        ca = ABS(qavg(i,5))  ! alfven speed
        cf = SQRT( MAX(0.D0, 0.5d0 * (tmp1 + tmp2) ))  ! fast mode speed
        cs = SQRT( MAX(0.D0, 0.5d0 * (tmp1 - tmp2) ))  ! slow mode speed
        
        ! wave speeds
        wave_speeds(i,1) = cs
        wave_speeds(i,2) = ca
        wave_speeds(i,3) = cf
        
        ! eigenvalues
        eigenvalues(i,1) = vx2 - cf
        eigenvalues(i,2) = vx2 - ca
        eigenvalues(i,3) = vx2 - cs
        eigenvalues(i,4) = vx2 + cs
        eigenvalues(i,5) = vx2 + ca
        eigenvalues(i,6) = vx2 + cf

    END DO 

    
    !PRINT*,'Wave Speeds: i, cs, ca, cf :'
    !DO i = 0, nx
    !   PRINT*,i, wave_speeds(i,1), wave_speeds(i,2), wave_speeds(i,3)
    !END DO
    
    
END SUBROUTINE compute_eigenvalues


! Compute eigenvectors
SUBROUTINE compute_eigenvectors()

    REAL*8 :: rho2, vx2, vy2, vz2, bx2, by2, bz2, a1, asqr, isrho, amc
    REAL*8 :: cf, ca, cs, bsqr, bt, tmp1, tmp2, tmp3, tmp4, tmp5, signbx, signbt, itheta1, itheta2
    REAL*8 :: alphas(1-nb:nx+nb), alphaf(1-nb:nx+nb), betay(1-nb:nx+nb), betaz(1-nb:nx+nb)
    REAL*8 :: delq(6), Lminus(6), Lplus(6)
    INTEGER :: i
      
    DO i = 1-nb, nx+nb-1

        asqr   = sound_speed**2      
        bt = SQRT(qavg(i,6)**2 + qavg(i,7)**2) ! magnitude of B-transverse
        
        ! compute eigenvector renormalization co-efficients, taking limiting values for degenracies (K98 eqn. 2.13-15)
        IF(bt .LT. ZERO .AND. ABS(asqr - wave_speeds(i,2)**2) .LT. ZERO) THEN
            alphaf(i) = 1.D0
            alphas(i) = 1.D0
        ELSE
            tmp1 = 1.D0 / SQRT(wave_speeds(i,3)**2 - wave_speeds(i,1)**2)
            alphaf(i) = SQRT(wave_speeds(i,3)**2 - wave_speeds(i,2)**2) * tmp1       
            alphas(i) = SQRT(wave_speeds(i,3)**2 - asqr) * tmp1
        END IF
 
        IF(bt .LT. ZERO) THEN
            betay(i) = 1.D0 / SQRT(2.D0)
            betaz(i) = 1.D0 / SQRT(2.D0)  
        ELSE
            tmp2 = 1.D0 / bt
            betay(i) = qavg(i,6) * tmp2
            betaz(i) = qavg(i,7) * tmp2       
        END IF
        
    END DO    
    
    !*************************************************
    ! compute the right eigenvectors (K98 eq 2.12)
    !*************************************************
    
    DO i = 1-nb, nx+nb-1
    
        isrho = 1.D0 / SQRT(qavg(i,1))       
        bt = SQRT(qavg(i,6)**2 + qavg(i,7)**2) ! magnitude of B-transverse
        amc = sound_speed**2 - wave_speeds(i,2)**2
        tmp1 = alphaf(i) * qavg(i,3)
        tmp2 = alphas(i) * betay(i) * qavg(i,5)
        tmp3 = alphaf(i) * qavg(i,4)
        tmp4 = alphas(i) * betaz(i) * qavg(i,5)
        tmp5 = alphas(i) * wave_speeds(i,3) * isrho
        
        ! fast mode minus eigenvector 
        Rk(i,1,1) = alphaf(i)
        Rk(i,2,1) = alphaf(i) * eigenvalues(i,1) 
        Rk(i,3,1) = tmp1 + tmp2
        Rk(i,4,1) = tmp3 + tmp4 
        Rk(i,5,1) = tmp5 * betay(i)
        Rk(i,6,1) = tmp5 * betaz(i)
        
        ! fast mode plus eigenvector 
        Rk(i,1,6) = alphaf(i)
        Rk(i,2,6) = alphaf(i) * eigenvalues(i,6) 
        Rk(i,3,6) = tmp1 - tmp2
        Rk(i,4,6) = tmp3 - tmp4 
        Rk(i,5,6) = tmp5 * betay(i)
        Rk(i,6,6) = tmp5 * betaz(i)
        

        ! check for continuity (K98 eqn. 2.20)
        IF(amc .LT. ZERO) THEN
        
            IF(qavg(i,6) .GE. ZERO .AND. qavg(i,7) .GT. ZERO) THEN
                signbt = 1.d0
            ELSE IF(qavg(i,6) .LE. ZERO .AND. qavg(i,7) .LT. ZERO) THEN
               signbt = -1.D0
            ELSE
               signbt = 1.d0            
            END IF
            
            Rk(i,:,1) = signbt * Rk(i,:,1) 
            Rk(i,:,6) = signbt * Rk(i,:,6)
        
        END IF


    END DO 

    DO i = 1-nb, nx+nb-1
    
        isrho = 1.D0 / SQRT(qavg(i,1))       
        signbx = SIGN(1.D0, qavg(i,5))
        
        ! alfven mode minus eigenvector 
        Rk(i,1,2) = 0.D0
        Rk(i,2,2) = 0.D0 
        Rk(i,3,2) = signbx * betaz(i)
        Rk(i,4,2) = -signbx * betay(i)
        Rk(i,5,2) = isrho * betaz(i)
        Rk(i,6,2) = -isrho * betay(i)
        
        ! alfven mode plus eigenvector 
        Rk(i,1,5) = 0.D0
        Rk(i,2,5) = 0.D0 
        Rk(i,3,5) = -signbx * betaz(i)
        Rk(i,4,5) = signbx * betay(i)
        Rk(i,5,5) = isrho * betaz(i)
        Rk(i,6,5) = -isrho * betay(i)
        
    END DO 
    
    DO i = 1-nb, nx+nb-1
    
        isrho = 1.D0 / SQRT(qavg(i,1))       
        signbx = SIGN(1.D0, qavg(i,5))
        amc = sound_speed**2 - wave_speeds(i,2)**2
        tmp1 = alphas(i) * qavg(i,3)
        tmp2 = alphaf(i) * betay(i) * sound_speed * signbx 
        tmp3 = alphas(i) * qavg(i,4)
        tmp4 = alphaf(i) * betaz(i) * sound_speed * signbx
        tmp5 = alphaf(i) * (sound_speed**2) * isrho /  wave_speeds(i,3)
        
        ! slow mode minus eigenvector 
        Rk(i,1,3) = alphas(i)
        Rk(i,2,3) = alphas(i) * eigenvalues(i,3) 
        Rk(i,3,3) = tmp1 - tmp2
        Rk(i,4,3) = tmp3 - tmp4 
        Rk(i,5,3) = -tmp5 * betay(i)
        Rk(i,6,3) = -tmp5 * betaz(i)
        
        ! slow mode plus eigenvector 
        Rk(i,1,4) = alphas(i)
        Rk(i,2,4) = alphas(i) * eigenvalues(i,4) 
        Rk(i,3,4) = tmp1 + tmp2
        Rk(i,4,4) = tmp3 + tmp4 
        Rk(i,5,4) = -tmp5 * betay(i)
        Rk(i,6,4) = -tmp5 * betaz(i)
        
        ! check for continuity (K98 eqn. 2.20)
        IF(amc .GT. ZERO) THEN
        
            IF(qavg(i,6) .GE. ZERO .AND. qavg(i,7) .GT. ZERO) THEN
                signbt = 1.d0
            ELSE IF(qavg(i,6) .LE. ZERO .AND. qavg(i,7) .LT. ZERO) THEN
               signbt = -1.D0
            ELSE
               signbt = 1.d0            
            END IF
        
            Rk(i,:,3) = signbt * Rk(i,:,3) 
            Rk(i,:,4) = signbt * Rk(i,:,4)
        
        END IF
        
    END DO 
    
    
    !*****************************************************
    ! compute the characteristic variables (K98 eqn. 2.27)
    !*****************************************************
    
    DO i = 1-nb, nx+nb-1
              
        signbx = SIGN(1.D0, qavg(i,5))
        amc = sound_speed**2 - wave_speeds(i,2)**2

        ! qavg_i+1 - qavg_i (Note: The Bx component has to be excluded)
        delq(1:4) = q(i+1,1:4) - q(i,1:4) 
        delq(5) = q(i+1,6) - q(i,6)
        delq(6) = q(i+1,7) - q(i,7)
 
        itheta1 = 1.D0 / (2.D0 * (alphaf(i)**2 * sound_speed**2 + alphas(i)**2 * wave_speeds(i,3)**2))  
        itheta2 = 1.D0 / (2.D0 * (alphaf(i)**2 * sound_speed * wave_speeds(i,3) + &
                  alphas(i)**2 * wave_speeds(i,1) * wave_speeds(i,2) ) )
       
        tmp1 = itheta1 * alphaf(i) * sound_speed**2
        tmp2 = itheta2 * ( -alphaf(i)*sound_speed*qavg(i,2) + &
               alphas(i)*wave_speeds(i,1)*signbx*(betay(i)*qavg(i,3)+betaz(i)*qavg(i,4)) )         
 
        tmp3 = itheta2 * alphas(i) * wave_speeds(i,1) * signbx 
        tmp4 = itheta1 * alphas(i) * wave_speeds(i,3) * SQRT(qavg(i,1))    
        
        ! Fast mode minus and plus eigenvector components (K98 eqn. 2.16-19)
        Lminus(1) =  tmp1 - tmp2
        Lplus(1)  = tmp1 + tmp2
        
        Lminus(2) = -itheta2 * alphaf(i) * sound_speed  
        Lplus(2)  = -Lminus(2)
       
        Lminus(3) = tmp3 * betay(i)   
        Lplus(3)  = -Lminus(3)
       
        Lminus(4) = tmp3 * betaz(i)
        Lplus(4)  = -Lminus(4)
       
        Lminus(5) = tmp4 * betay(i)
        Lplus(5)  = Lminus(5)     
        
        Lminus(6) = tmp4 * betaz(i)
        Lplus(6)  = Lminus(6)     
        
         ! check for continuity (K98 eqn. 2.20)
        IF(amc .LT. ZERO) THEN
        
            IF(qavg(i,6) .GE. ZERO .AND. qavg(i,7) .GT. ZERO) THEN
                signbt = 1.d0
            ELSE IF(qavg(i,6) .LE. ZERO .AND. qavg(i,7) .LT. ZERO) THEN
               signbt = -1.D0
            ELSE
               signbt = 1.d0            
            END IF
        
            Lminus = signbt * Lminus
            Lplus  = signbt * Lplus
        
        END IF
        
        
        ! Fast mode minus and plus chacteristic variables
        ck(i,1) = delq(1) * Lminus(1) + delq(2) * Lminus(2) + delq(3) * Lminus(3) + & 
                  delq(4) * Lminus(4) + delq(5) * Lminus(5) + delq(6) * Lminus(6)  

        ck(i,6) = delq(1) * Lplus(1) + delq(2) * Lplus(2) + delq(3) * Lplus(3) + & 
                  delq(4) * Lplus(4) + delq(5) * Lplus(5) + delq(6) * Lplus(6)
        
    END DO 
    
    
    DO i = 1-nb, nx+nb-1
              
        signbx = SIGN(1.D0, qavg(i,5))

        ! qavg_i+1 - qavg_i (Note: The Bx component has to be excluded)
        delq(1:4) = q(i+1,1:4) - q(i,1:4) 
        delq(5) = q(i+1,6) - q(i,6)
        delq(6) = q(i+1,7) - q(i,7)
        
        tmp1 = 0.5d0 * (betaz(i) * qavg(i,3) - betay(i) * qavg(i,4)) * signbx
        tmp2 = 0.5d0 * signbx
        tmp3 = 0.5d0 * SQRT(qavg(i,1))
        
        ! Alfven mode minus and plus eigenvector components (K98 eqn. 2.16-19)
        Lminus(1) =  -tmp1
        Lplus(1)  = tmp1 
        
        Lminus(2) = 0.D0  
        Lplus(2)  = 0.D0
       
        Lminus(3) = tmp2 * betaz(i)
        Lplus(3)  = -Lminus(3)
       
        Lminus(4) = -tmp2 * betay(i)
        Lplus(4)  = -Lminus(4)
       
        Lminus(5) = tmp3 * betaz(i)
        Lplus(5)  = Lminus(5)     
        
        Lminus(6) = -tmp3 * betay(i)
        Lplus(6)  = Lminus(6)     
        
        
        ! Alfven mode minus and plus chacteristic variables       
        ck(i,2) = delq(1) * Lminus(1) + delq(2) * Lminus(2) + delq(3) * Lminus(3) + & 
                  delq(4) * Lminus(4) + delq(5) * Lminus(5) + delq(6) * Lminus(6)  

        ck(i,5) = delq(1) * Lplus(1) + delq(2) * Lplus(2) + delq(3) * Lplus(3) + & 
                  delq(4) * Lplus(4) + delq(5) * Lplus(5) + delq(6) * Lplus(6)
                  
    END DO 
    
    
    DO i = 1-nb, nx+nb-1
              
        signbx = SIGN(1.D0, qavg(i,5))
        amc = sound_speed**2 - wave_speeds(i,2)**2

        ! qavg_i+1 - qavg_i (Note: The Bx component has to be excluded)
        delq(1:4) = q(i+1,1:4) - q(i,1:4) 
        delq(5) = q(i+1,6) - q(i,6)
        delq(6) = q(i+1,7) - q(i,7)
        
        itheta1 = 1.D0 / (2.D0 * (alphaf(i)**2 * sound_speed**2 + alphas(i)**2 * wave_speeds(i,3)**2))  
        itheta2 = 1.D0 / (2.D0 * (alphaf(i)**2 * sound_speed * wave_speeds(i,3) + &
                  alphas(i)**2 * wave_speeds(i,1) * wave_speeds(i,2) ) )
       
        tmp1 = itheta1 * alphas(i) * wave_speeds(i,3)**2
        tmp2 = itheta2 * ( alphas(i)*wave_speeds(i,2)*qavg(i,2) + &
               alphaf(i)*wave_speeds(i,3)*signbx*(betay(i)*qavg(i,3)+betaz(i)*qavg(i,4)) )         
 
        tmp3 = itheta2 * alphaf(i) * wave_speeds(i,3) * signbx 
        tmp4 = itheta1 * alphaf(i) * wave_speeds(i,3) * SQRT(qavg(i,1))    
        
        ! Slow mode minus and plus eigenvector components (K98 eqn. 2.16-19)
        Lminus(1) =  tmp1 + tmp2
        Lplus(1)  = tmp1 - tmp2
        
        Lminus(2) = -itheta2 * alphas(i) * wave_speeds(i,2)  
        Lplus(2)  = -Lminus(2)
       
        Lminus(3) = -tmp3 * betay(i)   
        Lplus(3)  = -Lminus(3)
       
        Lminus(4) = -tmp3 * betaz(i)
        Lplus(4)  = -Lminus(4)
       
        Lminus(5) = -tmp4 * betay(i)
        Lplus(5)  = Lminus(5)     
        
        Lminus(6) = -tmp4 * betaz(i)
        Lplus(6)  = Lminus(6)     
        
         ! check for continuity (K98 eqn. 2.20)
        IF(amc .GT. ZERO) THEN
        
            IF(qavg(i,6) .GE. ZERO .AND. qavg(i,7) .GT. ZERO) THEN
                signbt = 1.d0
            ELSE IF(qavg(i,6) .LE. ZERO .AND. qavg(i,7) .LT. ZERO) THEN
               signbt = -1.D0
            ELSE
               signbt = 1.d0            
            END IF
            
            Lminus = signbt * Lminus
            Lplus  = signbt * Lplus
        
        END IF
        
        
        ! Slow mode minus and plus chacteristic variables
        ck(i,3) = delq(1) * Lminus(1) + delq(2) * Lminus(2) + delq(3) * Lminus(3) + & 
                  delq(4) * Lminus(4) + delq(5) * Lminus(5) + delq(6) * Lminus(6)  

        ck(i,4) = delq(1) * Lplus(1) + delq(2) * Lplus(2) + delq(3) * Lplus(3) + & 
                  delq(4) * Lplus(4) + delq(5) * Lplus(5) + delq(6) * Lplus(6)
                  
    END DO 
    
END SUBROUTINE compute_eigenvectors


! Compute the time-averaged cell-interface fluxes
SUBROUTINE compute_TVD_fluxes()

    INTEGER :: i, k
    REAL*8 :: dxdt, gammak, betak
    REAL*8 :: Qk, chi
    REAL*8 :: gtilde1, gtilde2, Qk1, Qk2, chi1, chi2, signg
    REAL*8 :: gk(1-nb:nx+nb, nwaves), F(0:nx+1, nwaves)
    REAL*8 :: rho, vx, vy, vz, Bx, By, Bz


    dxdt = dx/dt

    ! clear tvd flux array 
    Fs = 0.d0
    
    ! compute cell center fluxes 
    DO i = 0, nx+1
    
        rho = q(i,1)
        vx = q(i,2) / rho
        vy = q(i,3) / rho
        vz = q(i,4) / rho
        Bx = q(i,5)
        By = q(i,6)
        Bz = q(i,7)
    
        F(i,1) = rho * vx
        F(i,2) = rho * (vx**2 + sound_speed**2) + 0.5d0 * (By**2 + Bz**2 - Bx**2)
        F(i,3) = rho * vx * vy - Bx * By
        F(i,4) = rho * vx * vz - Bx * Bz
        F(i,5) = vx * By - vy * Bx
        F(i,6) = vx * Bz - vz * Bx   
    
    END DO

    ! compute the flux limiters
     DO k = 1, nwaves
        DO i = 0, nx+1
                
            chi1 = dtdx*eigenvalues(i-1,k) 
            chi2 = dtdx*eigenvalues(i,k) 
        
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
            gtilde1 = 0.5d0 * (Qk1 - chi1**2) * ck(i-1,k) 
            
            !gk_tilde_i+1/2 (K98 eqn. 2.30)
            gtilde2 = 0.5d0 * (Qk2 - chi2**2) * ck(i,k) 
         
            signg = SIGN(1.D0, gtilde2)             
            
            ! flux limiter (minmod limiter)
            gk(i,k) = signg * MAX( 0.d0, MIN(ABS(gtilde2), gtilde1*signg) )         
            
            !gk(i,k) = signg * MAX( 0.d0, MIN(ABS(gtilde2), limiter_beta*gtilde1*signg), MIN(limiter_beta*ABS(gtilde2), gtilde1*signg) )   ! Sweby limiter      
                
        END DO
    END DO
    

    DO k = 1, nwaves
        DO i = 0, nx
        
        
            IF(ABS(ck(i,k)) .GT. ZERO) THEN
                gammak = (gk(i+1,k) - gk(i,k)) / ck(i,k)
            ELSE
                gammak = 0.d0
            END IF
            
            chi = dtdx*eigenvalues(i,k) + gammak
       
            ! Qk_i+1/2 (K98 eqn. 2.31)
            IF(ABS(chi) .LT. (2.D0*eps(k))) THEN
                Qk = 0.25D0 * (chi**2 / eps(k)) + eps(k)
            ELSE
                Qk = ABS(chi)
            END IF   
            
            ! betak_i+1/2 (K98 eqn. 2.26)
            betak = Qk * ck(i,k) - (gk(i,k) + gk(i+1,k))  
            
            ! accumulate TVD flux contribution from each wave  (K98 eqn. 2.25)
            Fs(i,1) = Fs(i,1) - betak * Rk(i,1,k) 
            Fs(i,2) = Fs(i,2) - betak * Rk(i,2,k) 
            Fs(i,3) = Fs(i,3) - betak * Rk(i,3,k) 
            Fs(i,4) = Fs(i,4) - betak * Rk(i,4,k) 
            Fs(i,5) = Fs(i,5) - betak * Rk(i,5,k) 
            Fs(i,6) = Fs(i,6) - betak * Rk(i,6,k) 
            
        END DO
    END DO

    ! add remaining term in the roe TVD flux
    DO i = 0, nx
        Fs(i,1) = 0.5D0 * ( dxdt * Fs(i,1) + F(i,1) + F(i+1,1) )
        Fs(i,2) = 0.5D0 * ( dxdt * Fs(i,2) + F(i,2) + F(i+1,2) )
        Fs(i,3) = 0.5D0 * ( dxdt * Fs(i,3) + F(i,3) + F(i+1,3) )
        Fs(i,4) = 0.5D0 * ( dxdt * Fs(i,4) + F(i,4) + F(i+1,4) )
        Fs(i,5) = 0.5D0 * ( dxdt * Fs(i,5) + F(i,5) + F(i+1,5) )
        Fs(i,6) = 0.5D0 * ( dxdt * Fs(i,6) + F(i,6) + F(i+1,6) )   
    END DO
    
    !DO i = 0, nx
    !    WRITE(12,*) i, Fs(i,1), Fs(i,2), Fs(i,3),Fs(i,4),Fs(i,5),Fs(i,6)   
    !END DO
    
    
END SUBROUTINE compute_TVD_fluxes


! final volume update of the state vector
SUBROUTINE update_state_vector()

    INTEGER :: i
    
    ! update the conserved variables (Note: Bx, which is q(i,5), remains unchanged)
    DO i = 1, nx
        q(i,1) = q(i,1) - dtdx * (Fs(i,1) - Fs(i-1,1))
        q(i,2) = q(i,2) - dtdx * (Fs(i,2) - Fs(i-1,2))
        q(i,3) = q(i,3) - dtdx * (Fs(i,3) - Fs(i-1,3))
        q(i,4) = q(i,4) - dtdx * (Fs(i,4) - Fs(i-1,4))
        q(i,6) = q(i,6) - dtdx * (Fs(i,5) - Fs(i-1,5))
        q(i,7) = q(i,7) - dtdx * (Fs(i,6) - Fs(i-1,6))
    END DO


    ! check for unphysical densities and apply protections accordingly
    DO i = 1, nx
        IF(q(i,1) .LT. ZERO) THEN
            PRINT*,'Negative density detected in Cell # ',i,', RHO= ',q(i,1)
            STOP
        END IF        
    END DO


END SUBROUTINE update_state_vector


! compute maximum wave speed (need this for calculating the time step size)
SUBROUTINE get_max_wavespeed(max_speed)

    REAL*8, INTENT(INOUT) :: max_speed
    INTEGER :: i
    REAL*8 :: vx
    
    max_speed = 0.d0

    ! find the maximum fast mode speed on the grid
    DO i = 0, nx
  
       vx = ABS(q(i,2)) / q(i,1)
       max_speed = MAX(max_speed, vx + wave_speeds(i,3))
    
    END DO

    !IF(max_speed .LT. ZERO) THEN
    !    PRINT*,'ERROR!! ZERO SPEED!!'
    !    STOP
    !END IF


END SUBROUTINE get_max_wavespeed



END MODULE Iso_TVD_Solver_mod