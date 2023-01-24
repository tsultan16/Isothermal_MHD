PROGRAM FFT_simple

! References: Numerical Recipes, Press et' al' 


IMPLICIT NONE


INTEGER, PARAMETER :: N = 8 
REAL(8), PARAMETER :: twopi = 8.d0*ATAN(1.d0)
REAL(8), ALLOCATABLE :: fin(:), Fout(:)


IF( MOD(LOG(N)/LOG(2),2) .NE. 0) THEN
    PRINT*,'ERROR. N needs to be a power of 2...'
    STOP
END IF


! Allocate arrays. Note: The size is 2*N because 
! they will include both real and imaginary components, stored sucessively.
ALLOCATE(fin(1:2*N), Fout(1:2*N))

fin = 0.d0
Fout = 0.d0

! Initialize input data
CALL RANDOM_NUMBER(fin)

Fout = fin

CALL fft(Fout, N)





CONTAINS


! 1D Daniel-Lanczos FFT algorithm 
SUBROUTINE fft(dat, nn, sgn)

    REAL(8), INTENT(INOUT) :: dat(:)   ! input array gets replaced by the inverse DFT
    INTEGER, INTENT(IN) :: nn          ! size of DFT
    INTEGER, INTENT(IN) :: sgn         ! DFT for sgn = 1 and Inverse DFT for sgn = -1
    REAL(8) :: tempr, tempi, theta, wi, wr, wpi, wpr, wtemp
    INTEGER :: i, j, n, m, mmax, istep
    
  
    !***************************************
    ! Sort input array in bit-reversed order  
    !***************************************    
    n = 2* nn
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
        IF((m .GT. 2) .AND. (j .GT. m)) THEN
            j = j - m 
            m = m / 2
        END IF
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
    
    
    IF(sgn .EQ. -1) dat = dat / DBLE(N)



END SUBROUTINE fft()




END PROGRAM FFT_simple