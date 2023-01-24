PROGRAM DFT_simple


IMPLICIT NONE


INTEGER, PARAMETER :: nx = 128
REAL(8), PARAMETER :: twopi = 8.d0*DATAN(1.d0)
REAL(8), ALLOCATABLE :: fx_re(:), fx_im(:), fk_re(:), fk_im(:), fk2(:), fk2_re(:), fk2_im(:)
REAL(8) :: dx, x, k, Lx, tmp
INTEGER :: i, j

tmp = DLOG(DBLE(nx))/DLOG(DBLE(2))


IF( ABS(MOD(tmp,1.d0)) .GE. 1.d-20) THEN
    PRINT*,'ERROR. N needs to be a power of 2...'
    STOP
END IF


! Allocate arrays.
ALLOCATE(fx_re(0:nx-1), fx_im(0:nx-1), fk_re(0:nx-1), fk_im(0:nx-1), fk2(0:2*nx-1), fk2_re(0:nx-1), fk2_im(0:nx-1))

fx_re = 0.d0
fx_im = 0.d0
fk_re = 0.d0
fk_im = 0.d0
fk2_re = 0.d0
fk2_im = 0.d0
fk2 = 0.d0

Lx = 1.d0
dx = Lx / DBLE(nx)


! Initialize input data
DO i = 0, nx-1
    x = i * dx
    
    !IF(x .GE. 0.25d0 .AND. x .LE. 0.75d0) THEN
    !    fx_re(i) = 1.d0
    !ELSE
    !    fx_re(i) = 0.d0
    !END IF
   
    fx_re(i) = COS(4.d0*twopi*x)  
 
    
END DO


CALL dft(fx_re, fx_im, fk_re, fk_im, 1)

!fk_re = 0.d0
!fk_im = 0.d0

!DO i = 0, nx-1

!    j = i- nx/2
!    IF(j .EQ. -5 .OR. j .EQ. 5) fk_re(i) = 1.d10 
    
!END DO


fx_re = 0.d0
fx_im = 0.d0

CALL dft(fk_re, fk_im, fx_re, fx_im, -1)


!DO i = 0, 2*nx-1, 2
!    fk2(i) = fx_re(i/2)   ! real part
!    fk2(i+1) = fx_im(i/2) ! followed by imaginary part
!END DO

!CALL fft(fk2, 1)

!DO i = 0, 2*nx-1, 2
!    fk2_re(i/2) = fk2(i)   ! real part
!    fk2_im(i/2) = fk2(i+1) ! followed by imaginary part
!END DO


! rotate fk array
DO i = 0, nx/2-1

    tmp = fk_re(i)
    fk_re(i) = fk_re(i+nx/2)
    fk_re(i+nx/2) = tmp
    
    tmp = fk_im(i)
    fk_im(i) = fk_im(i+nx/2)
    fk_im(i+nx/2) = tmp
    
    !tmp = fk2_re(i)
    !fk2_re(i) = fk2_re(i+nx/2)
    !fk2_re(i+nx/2) = tmp
    
    !tmp = fk2_im(i)
    !fk2_im(i) = fk2_im(i+nx/2)
    !fk2_im(i+nx/2) = tmp
        
END DO



OPEN(UNIT=12, FILE='output.dat', FORM = 'UNFORMATTED', ACCESS = 'STREAM')

DO i = 0, nx-1
    k = i-nx/2 !* twopi / Lx
    x = i * dx
    WRITE(12) x, fx_re(i),k,fk_re(i),fk_im(i)!,fk2_re(i),fk2_im(i) 
END DO

CLOSE(UNIT=12)

DEALLOCATE(fx_re, fx_im, fk_re, fk_im, fk2_re, fk2_im, fk2)


PRINT*,'Done.'


CONTAINS


! This subroutine computes the 1D inverse DFT of an input sequence by direct summation  
SUBROUTINE dft(in_re, in_im, out_re, out_im, sgn)

    REAL(8), INTENT(IN) :: in_re(0:nx-1), in_im(0:nx-1)    
    REAL(8), INTENT(INOUT) :: out_re(0:nx-1), out_im(0:nx-1)   
    INTEGER, INTENT(IN) :: sgn                         ! DFT for sgn = 1 and Inverse DFT for sgn = -1
    INTEGER :: i, j
    REAL(8) :: theta, theta0
    
    
    theta0 =  twopi / DBLE(nx)
    
    ! clear output arrays
    out_re = 0.d0
    out_im = 0.d0
    
    DO i = 0, nx-1
        DO j = 0, nx-1
            theta = theta0 * i * j
            out_re(i) = out_re(i) + (in_re(j) * DCOS(theta) + sgn * in_im(j) * DSIN(theta)) 
            out_im(i) = out_im(i) + (in_im(j) * DCOS(theta) - sgn * in_re(j) * DSIN(theta))
        END DO
    END DO   
  
   IF(sgn .EQ. -1) THEN
       out_re = out_re / DBLE(nx)
       out_im = out_im / DBLE(nx)
   END IF


END SUBROUTINE dft




! 1D Daniel-Lanczos FFT algorithm 
SUBROUTINE fft(dat, sgn)

    REAL(8), INTENT(INOUT) :: dat(:)   ! input array gets replaced by output
    INTEGER, INTENT(IN) :: sgn         ! DFT for sgn = 1 and Inverse DFT for sgn = -1
    REAL(8) :: tempr, tempi, theta, wi, wr, wpi, wpr, wtemp
    INTEGER :: i, j, n, m, mmax, istep
    
  
    !***************************************
    ! Sort input array in bit-reversed order  
    !***************************************    
    n = 2* nx
    j = 1
    
    DO i = 1, n, 2
    
        IF(j .GT. i) THEN  ! swap the two complex numbers
            tempr = dat(j)
            tempi = dat(j+1)
            dat(j) = dat(i)
            dat(j+1) = dat(i+1)
            dat(i) = tempr
            dat(i+1) = tempi                    
        END IF
    
        m = n/2
        DO WHILE ((m .GT. 2) .AND. (j .GT. m))
            j = j - m 
            m = m / 2
        END DO
        j = j + m
    END DO
      
      
    !********************************************************************************
    ! Using Danielson-Laczos lemma, compute the DFT by summing up the 1-pt base DFT's
    !********************************************************************************     
    mmax = 2
    
    DO WHILE(n .GT. mmax) 
    
        ! initialize for trigonometric recurrence
        istep = 2 * mmax
        theta = twopi / (sgn*mmax)  
        wpr = -2.d0 * SIN(0.5d0 * theta)**2 
        wpi =  SIN(theta)
        wr = 1.d0
        wi = 0.d0
        
        DO m = 1, mmax, 2 
            DO i = m, n, istep
       
                j = i + mmax
                
                ! apply Danielson-Lanczos lemma
                tempr = wr*dat(j) - wi*dat(j+1)
                tempi = wr*dat(j+1) + wi*dat(j)
                dat(j) = dat(i) - tempr 
                dat(j+1) = dat(i+1) - tempi 
                dat(i) = dat(i) + tempr 
                dat(i+1) = dat(i+1) + tempi 
                
            END DO
            wtemp = wr
            wr = wr*wpr - wi*wpi + wr
            wi = wi*wpr + wtemp*wpi + wi
            
        END DO
        mmax = istep
        
    END DO
    
    
    IF(sgn .EQ. -1) dat = dat / DBLE(nx)



END SUBROUTINE fft


END PROGRAM DFT_simple