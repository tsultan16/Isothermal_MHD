PROGRAM analysis

USE constants_mod
USE fft_mod
USE readfile_mod
IMPLICIT NONE

INCLUDE 'mpif.h'

!##########################################################################################

INTEGER, PARAMETER :: ndumps = 10
INTEGER, PARAMETER :: nbegin = 0 

REAL*4, ALLOCATABLE :: fx(:,:,:,:), fft_in(:,:,:), vp(:,:,:,:)
REAL*4, ALLOCATABLE :: vxk_re(:,:,:), vxk_im(:,:,:),vyk_re(:,:,:), vyk_im(:,:,:),vzk_re(:,:,:), vzk_im(:,:,:), &
                       bxk_re(:,:,:), bxk_im(:,:,:),byk_re(:,:,:), byk_im(:,:,:),bzk_re(:,:,:), bzk_im(:,:,:), &
                       rhok_re(:,:,:), rhok_im(:,:,:), va_k_re(:,:,:), va_k_im(:,:,:), vf_k_re(:,:,:), vf_k_im(:,:,:), &
                       vs_k_re(:,:,:), vs_k_im(:,:,:)
                       
REAL*4 :: Lx, Ly, Lz, dx
REAL*4 :: x, y, z, rho, bmag, bhat(3), vdotb
REAL*4 :: total_ke, v_rms, b_rms
INTEGER :: nmax, i, j, k, nd, mem_bytes
REAL*4 :: t0, t1, t2, t2a, t3, t4, t5, t6, t7, t8, t9, t10

!##########################################################################################
nmax = MIN(nranks_x*nx,nranks_y*ny,nranks_z*nz)

ALLOCATE(vxk_re(0:nmax-1,0:nmax-1,0:nmax-1))
ALLOCATE(vxk_im(0:nmax-1,0:nmax-1,0:nmax-1))
ALLOCATE(vyk_re(0:nmax-1,0:nmax-1,0:nmax-1))
ALLOCATE(vyk_im(0:nmax-1,0:nmax-1,0:nmax-1))
ALLOCATE(vzk_re(0:nmax-1,0:nmax-1,0:nmax-1))
ALLOCATE(vzk_im(0:nmax-1,0:nmax-1,0:nmax-1))
ALLOCATE(bxk_re(0:nmax-1,0:nmax-1,0:nmax-1))
ALLOCATE(bxk_im(0:nmax-1,0:nmax-1,0:nmax-1))
ALLOCATE(byk_re(0:nmax-1,0:nmax-1,0:nmax-1))
ALLOCATE(byk_im(0:nmax-1,0:nmax-1,0:nmax-1))
ALLOCATE(bzk_re(0:nmax-1,0:nmax-1,0:nmax-1))
ALLOCATE(bzk_im(0:nmax-1,0:nmax-1,0:nmax-1))
ALLOCATE(rhok_re(0:nmax-1,0:nmax-1,0:nmax-1))
ALLOCATE(rhok_im(0:nmax-1,0:nmax-1,0:nmax-1))
ALLOCATE(va_k_re(0:nmax-1,0:nmax-1,0:nmax-1))
ALLOCATE(va_k_im(0:nmax-1,0:nmax-1,0:nmax-1))
ALLOCATE(vf_k_re(0:nmax-1,0:nmax-1,0:nmax-1))
ALLOCATE(vf_k_im(0:nmax-1,0:nmax-1,0:nmax-1))
ALLOCATE(vs_k_re(0:nmax-1,0:nmax-1,0:nmax-1))
ALLOCATE(vs_k_im(0:nmax-1,0:nmax-1,0:nmax-1))

ALLOCATE(fx(1:nx*nranks_x,1:ny*nranks_y,1:nz*nranks_z,7))
ALLOCATE(vp(1:nx*nranks_x,1:ny*nranks_y,1:nz*nranks_z,3))
ALLOCATE(fft_in(1:nx*nranks_x,1:ny*nranks_y,1:nz*nranks_z))

mem_bytes = 20*SIZEOF(vxk_re) + SIZEOF(fx) + SIZEOF(fft_in) + SIZEOF(vp)
 
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
DO nd = nbegin, ndumps

    t1 = MPI_WTIME()
    
    fx = 0.d0
    
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

    t2a = MPI_WTIME()


    PRINT*,'Performing FFT of density.'
    
    ! perform FFTs on velocity field
    PRINT*,'FFT(vx)'
    fft_in(:,:,:) = fx(:,:,:,1) 
    CALL fft_3d(fft_in, vxk_re, vxk_im)
    vyk_re = 0.d0
    vyk_im = 0.d0
    vzk_re = 0.d0
    vzk_im = 0.d0
    
    PRINT*,'Computing density power spectrum.'
    CALL compute_pk(nd,0,vxk_re,vxk_im,vyk_re,vyk_im,vzk_re,vzk_im)



    PRINT*,'Performing FFT of velocity field.'

    t3 = MPI_WTIME()

    ! perform FFTs on velocity field
    PRINT*,'FFT(vx)'
    fft_in(:,:,:) = fx(:,:,:,2) 
    CALL fft_3d(fft_in, vxk_re, vxk_im)
    PRINT*,'FFT(vy)'
    fft_in(:,:,:) =fx(:,:,:,3) 
    CALL fft_3d(fft_in, vyk_re, vyk_im)
    PRINT*,'FFT(vz)'
    fft_in(:,:,:) = fx(:,:,:,4) 
    CALL fft_3d(fft_in, vzk_re, vzk_im)

    PRINT*,'Computing velocity power spectrum.'
    CALL compute_pk(nd,1,vxk_re,vxk_im,vyk_re,vyk_im,vzk_re,vzk_im)

    !CALL shift_fft(fxk_re, fxk_im)
    !CALL save_fft_to_file(nd)

    t4 = MPI_WTIME()

    
    PRINT*,'Performing FFT of magnetic field.'

    ! perform FFTs on magnetic field
    PRINT*,'FFT(bx)'
    fft_in(:,:,:) = fx(:,:,:,5) 
    CALL fft_3d(fft_in, bxk_re, bxk_im)
    PRINT*,'FFT(by)'
    fft_in(:,:,:) =fx(:,:,:,6) 
    CALL fft_3d(fft_in, byk_re, byk_im)
    PRINT*,'FFT(bz)'
    fft_in(:,:,:) = fx(:,:,:,7) 
    CALL fft_3d(fft_in, bzk_re, bzk_im)

    PRINT*,'Computing Bfield power spectrum.'
    CALL compute_pk(nd,2,bxk_re,bxk_im,byk_re,byk_im,bzk_re,bzk_im)

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

    ! perform FFTs on parallel velocity field
    PRINT*,'FFT(vparx)'
    fft_in(:,:,:) = vp(:,:,:,1) 
    CALL fft_3d(fft_in, vxk_re, vxk_im)
    PRINT*,'FFT(vpary)'
    fft_in(:,:,:) = vp(:,:,:,2) 
    CALL fft_3d(fft_in, vyk_re, vyk_im)
    PRINT*,'FFT(vparz)'
    fft_in(:,:,:) = vp(:,:,:,3) 
    CALL fft_3d(fft_in, vzk_re, vzk_im)

    PRINT*,'Computing v parallel power spectrum.'
    CALL compute_pk(nd,3,vxk_re,vxk_im,vyk_re,vyk_im,vzk_re,vzk_im)
    
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

    ! perform FFTs on perp. velocity field
    PRINT*,'FFT(vperpx)'
    fft_in(:,:,:) = vp(:,:,:,1) 
    CALL fft_3d(fft_in, vxk_re, vxk_im)
    PRINT*,'FFT(vperpy)'
    fft_in(:,:,:) = vp(:,:,:,2) 
    CALL fft_3d(fft_in, vyk_re, vyk_im)
    PRINT*,'FFT(vperpz)'
    fft_in(:,:,:) = vp(:,:,:,3) 
    CALL fft_3d(fft_in, vzk_re, vzk_im)

    PRINT*,'Computing v perp power spectrum.'
    CALL compute_pk(nd,4,vxk_re,vxk_im,vyk_re,vyk_im,vzk_re,vzk_im)


    t9 = MPI_WTIME()
    !##############################################################################################################
    PRINT*,'Performing MHD mode decomposition..'
    
    CALL mode_decomposition(nd)

    t10 = MPI_WTIME()
    !##############################################################################################################


    PRINT*,''
    PRINT*,'File read time (sec) = ', t2-t1
    PRINT*,'FFT time (sec) = ', t3-t2a+t5-t3+t7-t6+t9-t8
    PRINT*,'Mode decomposition time (sec)= ', t10-t9
    WRITE(*,'(" Time elapsed = ",i3," hour ",i3," min ",i3, "seconds.")') INT(t9-t0)/3600 , MOD(INT(t9-t0),3600)/60, MOD(INT(t9-t0),60) 
    PRINT*,''

END DO





CLOSE(UNIT=3)

DEALLOCATE(fx, vp, fft_in)
DEALLOCATE(vxk_re, vxk_im, vyk_re, vyk_im, vzk_re, vzk_im)
DEALLOCATE(bxk_re, bxk_im, byk_re, byk_im, bzk_re, bzk_im)
DEALLOCATE(rhok_re, rhok_im)
DEALLOCATE(va_k_re, va_k_im, vf_k_re, vf_k_im, vs_k_re, vs_k_im)

PRINT*,'Done.'
PRINT*,''


!##########################################################################################


CONTAINS


! fourier MHD wave mode decomposition (Ref: Cho & Lazarian, MNRAS, 345, 2003)
SUBROUTINE mode_decomposition(t_dump)

    INTEGER, INTENT(IN) :: t_dump 
    REAL*4 :: bmag, B0(3), alpha, D, zet_a(3), zet_f(3), zet_s(3), khat(3), &
              kdotB, kvec(3), kpar(3), kperp(3), kpar_hat(3), kperp_hat(3), tmp1, tmp2, rhok_f(2), rhok_s(2), bk_a(2), bk_f(2), bk_s(2), &
              delvk_f, delvk_s, cf, cs, vsqr_a, vsqr_f, vsqr_s, rhosqr_f, rhosqr_s, bsqr_a, bsqr_f, bsqr_s  
    INTEGER :: kx, ky, kz, kbins, ik
    REAL*4 :: kmin, kmax, dk, k, k0, dv_shell 
    CHARACTER(LEN=300) :: filename
    CHARACTER(LEN=6) :: uniti
    REAL*4, ALLOCATABLE :: Pk_v(:,:), Pk_b(:,:), Pk_rho(:,:) 


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
    
    ! set number of bins (i.e. k-shells)
    kbins = nmax
    
    ALLOCATE(Pk_v(1:kbins,3),Pk_b(1:kbins,3),Pk_rho(1:kbins,2))
    Pk_v = 0.d0
    Pk_b = 0.d0
    Pk_rho = 0.d0
    
    ! k band parameters
    k0 = TWOPI/Lx
    kmin = 0.d0
    kmax = SQRT(3.d0)*k0*(nmax/2)    
    dk = (kmax-kmin)/(kbins) ! uinform bin width

    B0(:) = 0.d0

    ! compute global mean magnetic field
    DO k = 1, nranks_z*nz
        DO j = 1, nranks_y*ny
            DO i = 1, nranks_x*nx

              B0(1) = B0(1) + fx(i,j,k,5)
              B0(2) = B0(2) + fx(i,j,k,6)
              B0(3) = B0(3) + fx(i,j,k,7)
            
            END DO
        END DO
    END DO
    B0(:) = B0(:) /(nmax**3) 
    
    bmag = SQRT(SUM(B0**2))
    alpha = 1.d0/SUM(B0**2)


    PRINT*,''
    PRINT*,'<B0> = ',B0
    PRINT*,''

    ! loop over k space grid
    DO kz = 0, nmax/2-1
        DO ky = 0, nmax/2-1
            DO kx = 0, nmax/2-1

                IF((kx .EQ. 0) .AND. (ky .EQ. 0)  .AND. (kz .EQ. 0)) CYCLE

                k = k0 * SQRT(SNGL(kx**2 + ky**2 + kz**2))
                
                kvec(1) = k0*kx
                kvec(2) = k0*ky
                kvec(3) = k0*kz
                
                khat(:) = kvec(:)/k 

                ! paralell and perpendicular components of wave vector     
                kpar_hat(:) = B0(:) / bmag
                kpar(:) = (kvec(1)*kpar_hat(1)+kvec(2)*kpar_hat(2)+kvec(3)*kpar_hat(3)) * kpar_hat(:)
                kperp(:) = kvec(:) - kpar(:)
                
                kperp_hat(:) = 0.d0
                IF(SUM(kperp**2) .GT. 0.d0) kperp_hat(:) = kperp(:) / SQRT(SUM(kperp**2))
                
                ! compute displacement unit vectors for all three wave modes
                D = MAX(0.d0, (1.d0+alpha)**2 - 4.d0*alpha*SUM(kpar**2)/(k**2))              
    
                zet_a(1) = kperp_hat(2)*kpar_hat(3) - kperp_hat(3)*kpar_hat(2)
                zet_a(2) = kperp_hat(3)*kpar_hat(1) - kperp_hat(1)*kpar_hat(3)
                zet_a(3) = kperp_hat(1)*kpar_hat(2) - kperp_hat(2)*kpar_hat(1)
                
                zet_f(:) = (-1.d0+alpha+SQRT(D))*kpar(:) + (1.d0+alpha+SQRT(D))*kperp(:) 
                IF(SUM(zet_f**2) .GT. 0) zet_f(:) = zet_f(:)/(SQRT(SUM(zet_f**2)))
                              
                zet_s(:) = (-1.d0+alpha-SQRT(D))*kpar(:) + (1.d0+alpha-SQRT(D))*kperp(:)  
                IF(SUM(zet_s**2) .GT. 0) zet_s(:) = zet_s(:)/(SQRT(SUM(zet_s**2)))
                
                ! project velocity fourier components onto displacement unit vectors of each mode 
                
                ! alfven mode
                va_k_re(kx,ky,kz) = vxk_re(kx,ky,kz)*zet_a(1) + vyk_re(kx,ky,kz)*zet_a(2) + vzk_re(kx,ky,kz)*zet_a(3)
                va_k_im(kx,ky,kz) = vxk_im(kx,ky,kz)*zet_a(1) + vyk_im(kx,ky,kz)*zet_a(2) + vzk_im(kx,ky,kz)*zet_a(3)      
                
                ! fast mode
                vf_k_re(kx,ky,kz) = vxk_re(kx,ky,kz)*zet_f(1) + vyk_re(kx,ky,kz)*zet_f(2) + vzk_re(kx,ky,kz)*zet_f(3)  
                vf_k_im(kx,ky,kz) = vxk_im(kx,ky,kz)*zet_f(1) + vyk_im(kx,ky,kz)*zet_f(2) + vzk_im(kx,ky,kz)*zet_f(3) 
                
                ! slow mode
                vs_k_re(kx,ky,kz) = vxk_re(kx,ky,kz)*zet_s(1) + vyk_re(kx,ky,kz)*zet_s(2) + vzk_re(kx,ky,kz)*zet_s(3)                 
                vs_k_im(kx,ky,kz) = vxk_im(kx,ky,kz)*zet_s(1) + vyk_im(kx,ky,kz)*zet_s(2) + vzk_im(kx,ky,kz)*zet_s(3)                 
               
            END DO
        END DO
    END DO

    ! loop over k space grid and deposit fourier amplitudes (squared) into power spectrum bins (i.e. shells in k-space)
     DO kz = 0, nmax/2-1
        DO ky = 0, nmax/2-1
            DO kx = 0, nmax/2-1
                
                IF((kx .EQ. 0) .AND. (ky .EQ. 0)  .AND. (kz .EQ. 0)) CYCLE

                k = k0 * SQRT(SNGL(kx**2 + ky**2 + kz**2))
                
                kvec(1) = k0*kx
                kvec(2) = k0*ky
                kvec(3) = k0*kz
                
                khat(:) = kvec(:)/k 

                ! paralell and perpendicular components of wave vector     
                kpar_hat(:) = B0(:) / bmag
                kpar(:) = (kvec(1)*kpar_hat(1)+kvec(2)*kpar_hat(2)+kvec(3)*kpar_hat(3)) * kpar_hat(:)
                kperp(:) = kvec(:) - kpar(:)
                
                kperp_hat(:) = 0.d0
                IF(SUM(kperp**2) .GT. 0.d0) kperp_hat(:) = kperp(:) / SQRT(SUM(kperp**2))
                
               ! compute displacement unit vectors for all three wave modes
                D = MAX(0.d0, (1.d0+alpha)**2 - 4.d0*alpha*SUM(kpar**2)/(k**2))              
    
                zet_a(1) = kperp_hat(2)*kpar_hat(3) - kperp_hat(3)*kpar_hat(2)
                zet_a(2) = kperp_hat(3)*kpar_hat(1) - kperp_hat(1)*kpar_hat(3)
                zet_a(3) = kperp_hat(1)*kpar_hat(2) - kperp_hat(2)*kpar_hat(1)
                
                zet_f(:) = (-1.d0+alpha+SQRT(D))*kpar(:) + (1.d0+alpha+SQRT(D))*kperp(:) 
                IF(SUM(zet_f**2) .GT. 0) zet_f(:) = zet_f(:)/(SQRT(SUM(zet_f**2)))
                              
                zet_s(:) = (-1.d0+alpha-SQRT(D))*kpar(:) + (1.d0+alpha-SQRT(D))*kperp(:)  
                IF(SUM(zet_s**2) .GT. 0) zet_s(:) = zet_s(:)/(SQRT(SUM(zet_s**2)))
                
                
                tmp1 = khat(1)*zet_f(1) + khat(2)*zet_f(2) + khat(3)*zet_f(3) 
                tmp2 = khat(1)*zet_s(1) + khat(2)*zet_s(2) + khat(3)*zet_s(3) 
                delvk_f = 2.d0 * vf_k_im(kx,ky,kz)
                delvk_s = 2.d0 * vs_k_im(kx,ky,kz)
                
                cf = SQRT(0.5d0 * (bmag**2) * (1.d0+alpha+SQRT(D)))
                cs = SQRT(MAX(0.d0, 0.5d0 * (bmag**2) * (1.d0+alpha-SQRT(D))))
                
                ! density fourier amplitudes (real part is zero)
                rhok_f(1) = 0.d0             
                rhok_s(1) = 0.d0
                rhok_f(2) = delvk_f * tmp1 /cf             
                rhok_s(2) = delvk_s * tmp2 /cs            
                
                
                tmp1 = SQRT((B0(2)*zet_f(3) - B0(3)*zet_f(2))**2 + (B0(3)*zet_f(1) - B0(1)*zet_f(3))**2 + &
                            (B0(1)*zet_f(2) - B0(2)*zet_f(1))**2)
                
                tmp2 = SQRT((B0(2)*zet_s(3) - B0(3)*zet_s(2))**2 + (B0(3)*zet_s(1) - B0(1)*zet_s(3))**2 + &
                            (B0(1)*zet_s(2) - B0(2)*zet_s(1))**2) 
                
                ! magnetic field fourier amplitudes
                bk_a(1) = va_k_re(kx,ky,kz)
                bk_a(2) = va_k_im(kx,ky,kz)
                bk_f(1) = 0.d0
                bk_f(2) = delvk_f * tmp1 /cf 
                bk_s(1) = 0.d0 
                bk_s(2) = delvk_s * tmp2 /cs 
            

                ! compute amplitude squared            
                vsqr_a = 2.d0 *(va_k_re(kx,ky,kz)**2 + va_k_im(kx,ky,kz)**2)
                vsqr_f = 2.d0 *(vf_k_re(kx,ky,kz)**2 + vf_k_im(kx,ky,kz)**2)
                vsqr_s = 2.d0 *(vs_k_re(kx,ky,kz)**2 + vs_k_im(kx,ky,kz)**2)
                
                rhosqr_f =  2.d0 *(rhok_f(1)**2 + rhok_f(2)**2)
                rhosqr_s =  2.d0 *(rhok_s(1)**2 + rhok_s(2)**2)
                
                
                bsqr_a =  2.d0 *(bk_a(1)**2 + bk_a(2)**2)
                bsqr_f =  2.d0 *(bk_f(1)**2 + bk_f(2)**2)
                bsqr_s =  2.d0 *(bk_s(1)**2 + bk_s(2)**2)
                
                ! deposit into k-bins
                ik = 1 + INT(k/dk)
                
                tmp1 = kmin + (ik-1) * dk 
                dv_shell = (FOURPI/3.d0)* ((tmp1+dk)**3 - tmp1**3)

                Pk_v(ik,1) = Pk_v(ik,1) + vsqr_a !* dv_shell
                Pk_v(ik,2) = Pk_v(ik,2) + vsqr_f !* dv_shell
                Pk_v(ik,3) = Pk_v(ik,3) + vsqr_s !* dv_shell
                Pk_b(ik,1) = Pk_b(ik,1) + bsqr_a !* dv_shell
                Pk_b(ik,2) = Pk_b(ik,2) + bsqr_f !* dv_shell
                Pk_b(ik,3) = Pk_b(ik,3) + bsqr_s !* dv_shell
                Pk_rho(ik,1) = Pk_rho(ik,1)  + rhosqr_f !* dv_shell
                Pk_rho(ik,2) = Pk_rho(ik,2)  + rhosqr_s !* dv_shell
    
            END DO
        END DO
    END DO    
        
    Pk_v = Pk_v * (dx**3)    
    Pk_b = Pk_b * (dx**3)    
    Pk_rho = Pk_rho * (dx**3)    
        
    ! dump power spectrum into file
    filename = TRIM('Output/modes_Pk_dump=')//TRIM(uniti)//TRIM('.dat')


    OPEN(UNIT=1, FILE=filename, FORM = 'UNFORMATTED', ACCESS = 'STREAM')
    
    DO ik = 1, kbins
        WRITE(1) (ik-1)*dk,Pk_v(ik,1),Pk_v(ik,2),Pk_v(ik,3),Pk_b(ik,1),Pk_b(ik,2),Pk_b(ik,3),Pk_rho(ik,1),Pk_rho(ik,2) 
    END DO
    
    CLOSE(UNIT=1)
    
    DEALLOCATE(Pk_v, Pk_b, Pk_rho)
 

END SUBROUTINE mode_decomposition


SUBROUTINE compute_pk(t_dump, fieldtype, fxk_re, fxk_im, fyk_re, fyk_im, fzk_re, fzk_im)

    INTEGER, INTENT(IN) :: t_dump 
    INTEGER, INTENT(IN) :: fieldtype 
    REAL*4, INTENT(IN) :: fxk_re(0:nmax-1,0:nmax-1,0:nmax-1), fxk_im(0:nmax-1,0:nmax-1,0:nmax-1), &
                          fyk_re(0:nmax-1,0:nmax-1,0:nmax-1), fyk_im(0:nmax-1,0:nmax-1,0:nmax-1), &
                          fzk_re(0:nmax-1,0:nmax-1,0:nmax-1), fzk_im(0:nmax-1,0:nmax-1,0:nmax-1)
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
    kbins = nmax
    
    ALLOCATE(Pk(1:kbins))
    Pk = 0.d0
    
    ! k band parameters
    k0 = TWOPI/Lx
    kmin = 0.d0
    kmax = SQRT(3.d0)*k0*(nmax/2)    
    dk = (kmax-kmin)/(kbins) ! uinform bin width
    
    ! loop over k space grid and deposit velocity fourier amplitudes (squared) into power spectrum bins (i.e. shells in k-space)
     DO kz = 0, nmax/2-1
        DO ky = 0, nmax/2-1
            DO kx = 0, nmax/2-1
                
                k = k0 * SQRT(SNGL(kx**2 + ky**2 + kz**2))
            
                fk_sqr = 2.d0 *(fxk_re(kx,ky,kz)**2 + fxk_im(kx,ky,kz)**2 + &
                         fyk_re(kx,ky,kz)**2 + fyk_im(kx,ky,kz)**2 + &
                         fzk_re(kx,ky,kz)**2 + fzk_im(kx,ky,kz)**2) 
                
                ! zero frequency component
                IF((kx .EQ. 0) .AND. (ky .EQ. 0) .AND. (kz .EQ. 0)) fk_sqr = 0.5d0 * fk_sqr
                
                ik = 1 + INT(k/dk)

                Pk(ik) = Pk(ik) + fk_sqr
    
            END DO
        END DO
    END DO
    
    ! add the zero frequency component separately
    !Pk(1) = fxk_re(0,0,0)**2 + fxk_im(0,0,0)**2 + fyk_re(0,0,0)**2 + fyk_im(0,0,0)**2 + fzk_re(0,0,0)**2 + fzk_im(0,0,0)**2 
    
    Pk = Pk * (dx**3)    
    
    ! dump power spectrum into file
    IF(fieldtype .EQ. 0) THEN
        filename = TRIM('Output/density_Pk_dump=')//TRIM(uniti)//TRIM('.dat')
    ELSE IF(fieldtype .EQ. 1) THEN
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
        WRITE(1) (ik-1)*dk,Pk(ik)
    END DO
    
    CLOSE(UNIT=1)
    
    DEALLOCATE(Pk)


END SUBROUTINE compute_pk



SUBROUTINE compute_density_pdf(t_dump)

    INTEGER, INTENT(IN) :: t_dump 
    INTEGER :: i, j, k, nbins, ibin
    REAL*4 :: rhomin, rhomax, drho, pdfnorm
    REAL*4, ALLOCATABLE :: PDF(:) 
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
    nbins = nmax
    
    ALLOCATE(PDF(1:nbins))
    PDF = 0.d0
    
    ! set upper/lower bounds and bin size
    rhomin = 0.d0
    rhomax = 4.d0    
    drho = (rhomax-rhomin)/(nbins) ! uinform bin width 
    
    
    DO k = 1, nranks_z*nz
        DO j = 1, nranks_y*ny
            DO i = 1, nranks_x*nx

                ibin = 1 + INT(fx(i,j,k,1)/drho)
                PDF(ibin) = PDF(ibin) + 1.d0
                
            END DO
        END DO
    END DO

    ! normalization
    pdfnorm = drho*SUM(PDF)
    PDF = PDF / pdfnorm

    ! dump density pdf into file
    filename = TRIM('Output/density_PDF_dump=')//TRIM(uniti)//TRIM('.dat')
    OPEN(UNIT=1, FILE=filename, FORM = 'UNFORMATTED', ACCESS = 'STREAM')
    
    DO ibin = 1, nbins
        WRITE(1) (ibin-0.5)*drho,PDF(ibin)
    END DO
    
    CLOSE(UNIT=1)
    
    DEALLOCATE(PDF)

END SUBROUTINE compute_density_pdf


!SUBROUTINE enstrophy()
!
!END SUBROUTINE enstrophy



SUBROUTINE shift_fft(fxk_re, fxk_im)

    REAL*4, INTENT(INOUT) :: fxk_re(-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2), &
                             fxk_im(-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2)  
    INTEGER :: kx, ky, kz, kx_shft, ky_shft, kz_shft
    REAL*4, ALLOCATABLE :: fs_re(:,:,:), fs_im(:,:,:)
    
    ALLOCATE(fs_re(-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2))
    ALLOCATE(fs_im(-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2,-nmax/2:-1+nmax/2))
    
    DO kz = -nmax/2, -1+nmax/2
        DO ky = -nmax/2, -1+nmax/2
            DO kx = -nmax/2, -1+nmax/2

                kx_shft = kx - nmax * (-1 + (SIGN(1,kx)))/2 
                ky_shft = ky - nmax * (-1 + (SIGN(1,ky)))/2 
                kz_shft = kz - nmax * (-1 + (SIGN(1,kz)))/2          

                fs_re(kx,ky,kz) = fxk_re(kx_shft,ky_shft,kz_shft)
                fs_im(kx,ky,kz) = fxk_im(kx_shft,ky_shft,kz_shft)
    
            END DO
        END DO
    END DO
    
    DO kz = -nmax/2, -1+nmax/2
        DO ky = -nmax/2, -1+nmax/2
            DO kx = -nmax/2, -1+nmax/2

                fxk_re(kx+nmax/2,ky+nmax/2,kz+nmax/2) = fs_re(kx,ky,kz)
                fxk_im(kx+nmax/2,ky+nmax/2,kz+nmax/2) = fs_im(kx,ky,kz)
    
            END DO
        END DO
    END DO
    
    

    DEALLOCATE(fs_re, fs_im)

END SUBROUTINE shift_fft



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



SUBROUTINE save_fft_to_file(t_dump)

    INTEGER, INTENT(IN) :: t_dump 
    INTEGER :: kx, ky, kz, kx_shft, ky_shft, kz_shft
    CHARACTER(LEN=300) :: filename
    CHARACTER(LEN=6) :: uniti


    PRINT*,'Saving FFT to file..'

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
    


    filename = TRIM('Output/velocity_vxk_dump=')//TRIM(uniti)//TRIM('.dat')

    OPEN(UNIT=1, FILE=filename, FORM = 'UNFORMATTED', ACCESS = 'STREAM')
    
    DO kz = -nmax/2, -1+nmax/2
        DO ky = -nmax/2, -1+nmax/2
            DO kx = -nmax/2, -1+nmax/2

                kx_shft = kx - nmax * (-1 + (SIGN(1,kx)))/2 
                ky_shft = ky - nmax * (-1 + (SIGN(1,ky)))/2 
                kz_shft = kz - nmax * (-1 + (SIGN(1,kz)))/2          

                !WRITE(1) fxk_re(kx_shft,ky_shft,kz_shft), fxk_im(kx_shft,ky_shft,kz_shft)
                !WRITE(1) fxk_re(kx+nmax/2,ky+nmax/2,kz+nmax/2), fxk_im(kx+nmax/2,ky+nmax/2,kz+nmax/2)
    
            END DO
        END DO
    END DO

    
    CLOSE(UNIT=1)

END SUBROUTINE save_fft_to_file


END PROGRAM analysis