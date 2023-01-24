PROGRAM analysis

USE constants_mod
USE fft_mod
USE readfile_mod
IMPLICIT NONE

INCLUDE 'mpif.h'

!##########################################################################################

INTEGER, PARAMETER :: ndumps = 4

REAL*4, ALLOCATABLE :: fx(:,:,:,:), fft_in(:,:,:), vp(:,:,:,:)
REAL*4, ALLOCATABLE :: fxk_re(:,:,:), fxk_im(:,:,:),fyk_re(:,:,:), fyk_im(:,:,:),fzk_re(:,:,:), fzk_im(:,:,:)
REAL*4 :: Lx, Ly, Lz, dx
REAL*4 :: x, y, z, rho, bmag, bhat(3), vdotb
REAL*4 :: total_ke, v_rms, b_rms
INTEGER :: nmax, i, j, k, nd, mem_bytes
REAL*4 :: t0, t1, t2, t3, t4, t5, t6, t7, t8, t9

!##########################################################################################
nmax = MIN(nranks_x*nx,nranks_y*ny,nranks_z*nz)

ALLOCATE(fxk_re(-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2))
ALLOCATE(fxk_im(-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2))
ALLOCATE(fyk_re(-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2))
ALLOCATE(fyk_im(-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2))
ALLOCATE(fzk_re(-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2))
ALLOCATE(fzk_im(-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2))
ALLOCATE(fx(1:nx*nranks_x,1:ny*nranks_y,1:nz*nranks_z,7))
ALLOCATE(vp(1:nx*nranks_x,1:ny*nranks_y,1:nz*nranks_z,3))
ALLOCATE(fft_in(1:nx*nranks_x,1:ny*nranks_y,1:nz*nranks_z))

mem_bytes = SIZEOF(fxk_re) + SIZEOF(fxk_im) + SIZEOF(fyk_re) + SIZEOF(fyk_im) + SIZEOF(fzk_re) + SIZEOF(fzk_im) +&
            SIZEOF(fx) + SIZEOF(fft_in) + SIZEOF(vp)
 
PRINT*,''
PRINT*,'Memory allocated for work arrays (Mb) = ', mem_bytes*1.e-6
PRINT*,''
!##########################################################################################

OPEN(UNIT=3, FILE='Output/avgs.dat', FORM = 'UNFORMATTED', ACCESS = 'STREAM')


! grid size
Lx = 1.d0
Ly = 1.d0
Lz = 1.d0
dx = Lx/DBLE(nranks_x*nx)


! initialize fft module
CALL fft_init(nranks_x*nx,nranks_y*ny,nranks_z*nz)

t0 = MPI_WTIME()

! loop over file dumps
DO nd = 0, ndumps


    t1 = MPI_WTIME()
    !CALL readfile_ranks_singlevar(0, 1, 8, .FALSE., fx)
    CALL readfile_ranks_multivar(nd, 1, 1, 7, .FALSE., fx)
    t2 =  MPI_WTIME()


    CALL compute_ke()
    CALL compute_vrms()
    CALL compute_brms()

    WRITE(3) SNGL(nd*dump_frequency),total_ke,v_rms, b_rms

    ! convert momentum into velocity
    PRINT*,'Converting momentum to velocity.'

    DO k = 1, nranks_z*nz
    DO j = 1, nranks_y*ny
    DO i = 1, nranks_x*nx

        rho = fx(i,j,k,1)
        fx(i,j,k,2) = fx(i,j,k,2)/rho  
        fx(i,j,k,3) = fx(i,j,k,3)/rho  
        fx(i,j,k,4) = fx(i,j,k,4)/rho  

    END DO
    END DO
    END DO


    PRINT*,'Performing FFT of velocity field.'

    t3 = MPI_WTIME()

    ! perform FFTs on velocity field
    PRINT*,'FFT(vx)'
    fft_in(:,:,:) = fx(:,:,:,2) 
    CALL fft_3d(fft_in, fxk_re, fxk_im)
    PRINT*,'FFT(vy)'
    fft_in(:,:,:) =fx(:,:,:,3) 
    CALL fft_3d(fft_in, fyk_re, fyk_im)
    PRINT*,'FFT(vz)'
    fft_in(:,:,:) = fx(:,:,:,4) 
    CALL fft_3d(fft_in, fzk_re, fzk_im)

    PRINT*,'Computing velocity power spectrum.'
    CALL compute_pk(nd,1)

    t4 = MPI_WTIME()

    GO TO 111
    
    
    PRINT*,'Performing FFT of magnetic field.'

    ! perform FFTs on velocity field
    PRINT*,'FFT(bx)'
    fft_in(:,:,:) = fx(:,:,:,5) 
    CALL fft_3d(fft_in, fxk_re, fxk_im)
    PRINT*,'FFT(by)'
    fft_in(:,:,:) =fx(:,:,:,6) 
    CALL fft_3d(fft_in, fyk_re, fyk_im)
    PRINT*,'FFT(bz)'
    fft_in(:,:,:) = fx(:,:,:,7) 
    CALL fft_3d(fft_in, fzk_re, fzk_im)

    PRINT*,'Computing Bfield power spectrum.'
    CALL compute_pk(nd,2)

    t5 = MPI_WTIME()

    !##############################################################################################################

    PRINT*,'Computing v parallel...'

    ! Compute velocity component parallel to local B field 
    DO k = 1, nranks_z*nz
    DO j = 1, nranks_y*ny
    DO i = 1, nranks_x*nx


        bmag = SQRT(fx(i,j,k,5)**2+fx(i,j,k,6)**2+fx(i,j,k,7)**2)
        bhat(1) = fx(i,j,k,5) /bmag
        bhat(2) = fx(i,j,k,6) /bmag
        bhat(3) = fx(i,j,k,7) /bmag
        vdotb = (fx(i,j,k,2)*fx(i,j,k,5) + fx(i,j,k,3)*fx(i,j,k,6) + fx(i,j,k,4)*fx(i,j,k,7))/bmag
        
        ! parallel component
        vp(i,j,k,1) = vdotb * bhat(1)   
        vp(i,j,k,2) = vdotb * bhat(2)   
        vp(i,j,k,3) = vdotb * bhat(3)   


    END DO
    END DO
    END DO
    
    t6 = MPI_WTIME()
     
    PRINT*,'Performing FFT of v parallel...'

    ! perform FFTs on velocity field
    PRINT*,'FFT(vparx)'
    fft_in(:,:,:) = vp(:,:,:,1) 
    CALL fft_3d(fft_in, fxk_re, fxk_im)
    PRINT*,'FFT(vpary)'
    fft_in(:,:,:) = vp(:,:,:,2) 
    CALL fft_3d(fft_in, fyk_re, fyk_im)
    PRINT*,'FFT(vparz)'
    fft_in(:,:,:) = vp(:,:,:,3) 
    CALL fft_3d(fft_in, fzk_re, fzk_im)

    PRINT*,'Computing v parallel power spectrum.'
    CALL compute_pk(nd,3)
    
    t7 = MPI_WTIME()

   !##############################################################################################################
    
    PRINT*,'Computing v perp...'

    ! Compute velocity component perpendicular to local B field 
    DO k = 1, nranks_z*nz
    DO j = 1, nranks_y*ny
    DO i = 1, nranks_x*nx


        bmag = SQRT(fx(i,j,k,5)**2+fx(i,j,k,6)**2+fx(i,j,k,7)**2)
        bhat(1) = fx(i,j,k,5) /bmag
        bhat(2) = fx(i,j,k,6) /bmag
        bhat(3) = fx(i,j,k,7) /bmag
        vdotb = (fx(i,j,k,2)*fx(i,j,k,5) + fx(i,j,k,3)*fx(i,j,k,6) + fx(i,j,k,4)*fx(i,j,k,7))/bmag 

        ! perpendicular component
        vp(i,j,k,1) = fx(i,j,k,2) - vdotb * bhat(1)   
        vp(i,j,k,2) = fx(i,j,k,3) - vdotb * bhat(2)    
        vp(i,j,k,3) = fx(i,j,k,4) - vdotb * bhat(3)   


    END DO
    END DO
    END DO

    t8 = MPI_WTIME()

    PRINT*,'Performing FFT of v perp...'

    ! perform FFTs on velocity field
    PRINT*,'FFT(vperpx)'
    fft_in(:,:,:) = vp(:,:,:,1) 
    CALL fft_3d(fft_in, fxk_re, fxk_im)
    PRINT*,'FFT(vperpy)'
    fft_in(:,:,:) = vp(:,:,:,2) 
    CALL fft_3d(fft_in, fyk_re, fyk_im)
    PRINT*,'FFT(vperpz)'
    fft_in(:,:,:) = vp(:,:,:,3) 
    CALL fft_3d(fft_in, fzk_re, fzk_im)

    PRINT*,'Computing v perp power spectrum.'
    CALL compute_pk(nd,4)

    111 CONTINUE

    t9 = MPI_WTIME()

   !##############################################################################################################


    PRINT*,''
    PRINT*,'File read time (sec) = ', t2-t1
    PRINT*,'FFT time (sec) = ', t5-t3+t7-t6+t9-t8
    WRITE(*,'(" Time elapsed = ",i3," hour ",i3," min ",i3, "seconds.")') INT(t9-t0)/3600 , MOD(INT(t9-t0),3600)/60, MOD(INT(t9-t0),60) 
    PRINT*,''

END DO





!OPEN(UNIT=12, FILE='output_3d.dat', FORM = 'UNFORMATTED', ACCESS = 'STREAM')

!DO iz = 1, nz
!DO iy = 1, ny
!DO ix = 1, nx
!    WRITE(12) fx(ix,iy,iz,1), fk_re(ix,iy,iz,1), fk_im(ix,iy,iz,1)
!END DO
!END DO
!END DO

CLOSE(UNIT=3)

DEALLOCATE(fx, vp, fxk_re, fxk_im, fyk_re, fyk_im, fzk_re, fzk_im, fft_in)




PRINT*,'Done.'
PRINT*,''


!##########################################################################################


CONTAINS


SUBROUTINE compute_pk(t_dump, fieldtype)

    INTEGER, INTENT(IN) :: t_dump 
    INTEGER, INTENT(IN) :: fieldtype 
    INTEGER :: kx, ky, kz, kbins, ik
    REAL*4 :: kmin, kmax, dk, k, k0, fk_sqr
    REAL*4, ALLOCATABLE :: Pk(:) 
    CHARACTER(LEN=300) :: filename
    CHARACTER(LEN=6) :: uniti


    IF(t_dump<10) THEN
        WRITE(uniti,'(I1.1)') t_dump
    ELSE IF(t_dump>=10 .and. t_dump<100) THEN
        WRITE(uniti,'(I2.2)') t_dump
    ELSE IF(t_dump>=100 .and. t_dump<1000) THEN
        WRITE (uniti,'(I3.3)') t_dump
    ELSE IF(t_dump>=1000 .and. t_dump<10000) THEN
        WRITE (uniti,'(I4.3)') t_dump
    ELSE IF(t_dump>=10000 .and. t_dump<100000) THEN
        WRITE (uniti,'(I5.3)') t_dump  
    END IF
    
    
    ! set number of bins 
    kbins = nmax/4
    
    ALLOCATE(Pk(1:kbins))
    Pk = 0.d0
    
    ! k band parameters
    k0 = TWOPI/Lx
    kmin = 0.d0
    kmax = SQRT(3.d0)*k0*(nmax/2)    
    dk = (kmax-kmin)/kbins ! uinform bin width
    
    ! loop over k space grid and deposit velocity fourier amplitudes (squared) into power spectrum bins (i.e. shells in k-space)
    DO kz = -(-1+nmax/2), -1+nmax/2
        DO ky = -(-1+nmax/2), -1+nmax/2
            DO kx = -(-1+nmax/2), -1+nmax/2
    
                k = k0 * SQRT(SNGL(kx**2 + ky**2 + kz**2))


                fk_sqr = fxk_re(kx,ky,kz)**2 + fxk_im(kx,ky,kz)**2 + fyk_re(kx,ky,kz)**2 + fyk_im(kx,ky,kz)**2 + fzk_re(kx,ky,kz)**2 + fzk_im(kx,ky,kz)**2 
                
                
                ik = 1 + INT(k/dk)

                Pk(ik) = Pk(ik) + fk_sqr
    
            END DO
        END DO
    END DO
    
        
    Pk = Pk * (dx**3)    
    
    ! dump power spectrum into file
    IF(fieldtype .EQ. 1) THEN
        filename = TRIM('Output/velocity_Pk_dump=')//TRIM(uniti)//TRIM('.dat')
    ELSE IF(fieldtype .EQ. 2) THEN
        filename = TRIM('Output/Bfield_Pk_dump=')//TRIM(uniti)//TRIM('.dat')    
    ELSE IF(fieldtype .EQ. 3) THEN
        filename = TRIM('Output/vpar_Pk_dump=')//TRIM(uniti)//TRIM('.dat')    
    ELSE IF(fieldtype .EQ. 4) THEN
        filename = TRIM('Output/vperp_Pk_dump=')//TRIM(uniti)//TRIM('.dat')    
    END IF    


    OPEN(UNIT=1, FILE=filename, FORM = 'UNFORMATTED', ACCESS = 'STREAM')
    
    DO ik = 1, kbins
        WRITE(1) (ik-0.5)*dk,Pk(ik)
    END DO
    
    CLOSE(UNIT=1)
    
    DEALLOCATE(Pk)


END SUBROUTINE compute_pk



SUBROUTINE compute_pk2(t_dump, fieldtype)

    INTEGER, INTENT(IN) :: t_dump 
    INTEGER, INTENT(IN) :: fieldtype 
    INTEGER :: kx, ky, kz, kbins, ik
    REAL*4 :: kmin, kmax, dk, k, k0, fk_sqr
    REAL*4, ALLOCATABLE :: Pk(:) 
    CHARACTER(LEN=300) :: filename
    CHARACTER(LEN=6) :: uniti


    IF(t_dump<10) THEN
        WRITE(uniti,'(I1.1)') t_dump
    ELSE IF(t_dump>=10 .and. t_dump<100) THEN
        WRITE(uniti,'(I2.2)') t_dump
    ELSE IF(t_dump>=100 .and. t_dump<1000) THEN
        WRITE (uniti,'(I3.3)') t_dump
    ELSE IF(t_dump>=1000 .and. t_dump<10000) THEN
        WRITE (uniti,'(I4.3)') t_dump
    ELSE IF(t_dump>=10000 .and. t_dump<100000) THEN
        WRITE (uniti,'(I5.3)') t_dump  
    END IF
    
    
    ! set number of bins 
    kbins = 2*nx*nranks_x
    
    ALLOCATE(Pk(1:kbins))
    Pk = 0.d0
    
    ! k band parameters
    k0 = TWOPI/Lx
    kmin = 0.d0
    kmax = 2.d0*k0*(nmax/2)    
    dk = (kmax-kmin)/kbins ! uinform bin width
    
    ! loop over k space grid and deposit velocity fourier amplitudes (squared) into power spectrum bins (i.e. shells in k-space)
    DO kz = 1, nmax
        DO ky = 1, nmax
            DO kx = 1, nmax
    
                k = k0 * SQRT(DBLE(kx**2 + ky**2 + kz**2))
                !fk_sqr = 2.d0 * (fxk_re(kx+nmax/2,ky+nmax/2,kz+nmax/2)**2 + fxk_im(kx+nmax/2,ky+nmax/2,kz+nmax/2)**2 + &
                !         fyk_re(kx+nmax/2,ky+nmax/2,kz+nmax/2)**2 + fyk_im(kx+nmax/2,ky+nmax/2,kz+nmax/2)**2 + &
                !         fzk_re(kx+nmax/2,ky+nmax/2,kz+nmax/2)**2 + fzk_im(kx+nmax/2,ky+nmax/2,kz+nmax/2)**2 )


                fk_sqr = fxk_re(kx+1+nmax/2,ky+1+nmax/2,kz+1+nmax/2)**2 + fxk_im(kx+1+nmax/2,ky+1+nmax/2,kz+1+nmax/2)**2 + &
                         fyk_re(kx+1+nmax/2,ky+1+nmax/2,kz+1+nmax/2)**2 + fyk_im(kx+1+nmax/2,ky+1+nmax/2,kz+1+nmax/2)**2 + &
                         fzk_re(kx+1+nmax/2,ky+1+nmax/2,kz+1+nmax/2)**2 + fzk_im(kx+1+nmax/2,ky+1+nmax/2,kz+1+nmax/2)**2 + &
                         fxk_re(-kx+1+nmax/2,-ky+1+nmax/2,-kz+1+nmax/2)**2 + fxk_im(-kx+1+nmax/2,-ky+1+nmax/2,-kz+1+nmax/2)**2 + &
                         fyk_re(-kx+1+nmax/2,-ky+1+nmax/2,-kz+1+nmax/2)**2 + fyk_im(-kx+1+nmax/2,-ky+1+nmax/2,-kz+1+nmax/2)**2 + &
                         fzk_re(-kx+1+nmax/2,-ky+1+nmax/2,-kz+1+nmax/2)**2 + fzk_im(-kx+1+nmax/2,-ky+1+nmax/2,-kz+1+nmax/2)**2
                         
                IF(kx .EQ. 0 .AND. ky .EQ. 0 .AND. kz .EQ. 0) fk_sqr = 0.5d0 * fk_sqr
    
                ik = 1 + INT(k/dk)

                Pk(ik) = Pk(ik) + fk_sqr
    
            END DO
        END DO
    END DO
    
    Pk = Pk * (dx**3)
    
    
    ! dump power spectrum into file
    IF(fieldtype .EQ. 1) THEN
        filename = TRIM('Output/velocity_Pk_dump=')//TRIM(uniti)//TRIM('.dat')
    ELSE IF(fieldtype .EQ. 2) THEN
        filename = TRIM('Output/Bfield_Pk_dump=')//TRIM(uniti)//TRIM('.dat')    
    ELSE IF(fieldtype .EQ. 3) THEN
        filename = TRIM('Output/vpar_Pk_dump=')//TRIM(uniti)//TRIM('.dat')    
    ELSE IF(fieldtype .EQ. 4) THEN
        filename = TRIM('Output/vperp_Pk_dump=')//TRIM(uniti)//TRIM('.dat')    
    END IF    


    OPEN(UNIT=1, FILE=filename, FORM = 'UNFORMATTED', ACCESS = 'STREAM')
    
    DO ik = 1, kbins
        WRITE(1) (ik-0.5)*dk,Pk(ik)
    END DO
    
    CLOSE(UNIT=1)
    
    DEALLOCATE(Pk)


END SUBROUTINE compute_pk2



SUBROUTINE compute_ke()

    INTEGER :: i, j, k
    REAL(4) :: rho 

    total_ke = 0.d0

    DO k = 1, nranks_z*nz
    DO j = 1, nranks_y*ny
    DO i = 1, nranks_x*nx

        rho = fx(i,j,k,1) 
        total_ke = total_ke + (fx(i,j,k,2)**2 + fx(i,j,k,3)**2 + fx(i,j,k,4)**2) /rho  

    END DO
    END DO
    END DO

    total_ke = 0.5d0 * total_ke * (dx**3)

    PRINT*,'Kinetic Energy = ',total_ke 


END SUBROUTINE compute_ke


SUBROUTINE compute_vrms()

    INTEGER :: i,j,k
    REAL(4) :: rho 

    v_rms = 0.d0

    DO k = 1, nranks_z*nz
    DO j = 1, nranks_y*ny
    DO i = 1, nranks_x*nx

        rho = fx(i,j,k,1) 
        v_rms = v_rms + (fx(i,j,k,2)**2 + fx(i,j,k,3)**2 + fx(i,j,k,4)**2) /(rho**2)  

    END DO
    END DO
    END DO

    v_rms = SQRT(v_rms*(dx**3))  

    PRINT*,'rms velocity = ',v_rms 


END SUBROUTINE compute_vrms


SUBROUTINE compute_brms()

    INTEGER :: i,j,k

    b_rms = 0.d0

    DO k = 1, nranks_z*nz
    DO j = 1, nranks_y*ny
    DO i = 1, nranks_x*nx

        b_rms = b_rms + (fx(i,j,k,5)**2 + fx(i,j,k,6)**2 + fx(i,j,k,7)**2)   

    END DO
    END DO
    END DO

    b_rms = SQRT(b_rms*(dx**3)) 

    PRINT*,'rms magnetic field = ',b_rms 


END SUBROUTINE compute_brms


END PROGRAM analysis