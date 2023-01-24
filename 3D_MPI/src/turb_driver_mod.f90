MODULE turb_driver_mod

USE constants_mod
USE grid_data_mod
USE mpi_domain_mod

IMPLICIT NONE

REAL*8, PARAMETER :: ONE = 1.d0 - 1.d-20
INTEGER :: quad_term1_req, quad_term2_req, momx_req, momy_req, momz_req
REAL*8 :: delv_rms, delv_rms1, delv_rms2, Lx, Ly, Lz, quad_term1, quad_term2, quad_term1_sum, quad_term2_sum, &
          momx_avg, momy_avg, momz_avg, momx_avg_sum, momy_avg_sum, momz_avg_sum 
REAL*8 :: khat(0:k_max,0:k_max,0:k_max,3)  
  
    
REAL*8, ALLOCATABLE :: rho0(:,:,:)   
    
INTEGER :: nnmax
    
CONTAINS


SUBROUTINE init_turb_driver(s1_in, s2_in, s3_in)

    INTEGER, INTENT(IN) :: s1_in, s2_in, s3_in
    INTEGER :: ix, iy, iz, kx, ky, kz, ixx, iyy, izz
    REAL*8 :: kx0, ky0, kz0, k
    INTEGER :: grid_memory

        
    ! assert if k_max value is legal
    IF(k_max .GE. MIN(nx*nranks_x, ny*nranks_y, nz*nranks_z)/2) THEN
        PRINT*,'ERROR!! Invalid k_max value.'
        STOP
    END IF
  
    nnmax = MAX(nx, ny, nz)
    
    ! set RNG seed values
    CALL set_taus_seed(s1_in, s2_in, s3_in)

    ! compute grid size
    Lx = DBLE(nx * nranks_x) * my_dx 
    Ly = DBLE(ny * nranks_y) * my_dx 
    Lz = DBLE(nz * nranks_z) * my_dx 

    ! precomupte k-space unit vectors
    kx0 = TWOPI / Lx
    ky0 = TWOPI / Ly
    kz0 = TWOPI / Lz 
    
    khat = 0.d0
    DO kz = 1, k_max
        DO ky = 1, k_max
            DO kx = 1, k_max 
                k = SQRT( (kx*kx0)**2 + (ky*ky0)**2  + (kz*kz0)**2 )
                khat(kx,ky,kz,1) = kx0 * kx / k
                khat(kx,ky,kz,2) = ky0 * ky / k
                khat(kx,ky,kz,3) = kz0 * kz / k
            END DO
        END DO
    END DO


    ! allocate memory for storing initial density
    ALLOCATE(rho0(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb))

    grid_memory =  SIZEOF(rho0) 
    
    
    IF(myrank .EQ. 0) THEN
        PRINT*,''
        PRINT*,'#####################################################################'
        PRINT*,'Total Memory for turbulence driver work array storage per MPI rank (Mb) = ', grid_memory*1.d-6                    
        PRINT*,'#####################################################################'
        PRINT*,''
    END IF
    
END SUBROUTINE init_turb_driver


SUBROUTINE destroy_turb_driver()

    DEALLOCATE(rho0)


END SUBROUTINE destroy_turb_driver


! set taus88 PRNG non-zero seed values
SUBROUTINE set_taus_seed(a, b, c)

    INTEGER, INTENT(IN) :: a, b, c

    s1 = a
    s2 = b
    s3 = c

END SUBROUTINE 



! returns a random number from [0,1)
FUNCTION rand_taus() RESULT(r)

    REAL*8:: r
    INTEGER :: b   ! Signed 32-bit INTEGER


    b = ISHFT(IEOR(ISHFT(s1,13),s1),-19)    
    s1 = IEOR(ISHFT(IAND(s1,4294967294),12),b)
    b = ISHFT(IEOR(ISHFT(s2,2),s2),-25)
    s2 = IEOR(ISHFT(IAND(s2,4294967288),4),b)
    b = ISHFT(IEOR(ISHFT(s3,3),s3),-11)
    s3 = IEOR(ISHFT(IAND(s3,4294967280),17),b)

    r = 0.5D0+IEOR(s1,IEOR(s2,s3))*2.3283064365D-10    ! 2.3283064365D-10 = 2^-32
    
    !r = MIN(r, ONE)
    
    !IF(r .GE. ONE .OR. r .LT. ZERO) THEN
    !    PRINT*,'Random number generator failed..r=',r
    !    STOP
    !END IF

END FUNCTION rand_taus



SUBROUTINE get_rho_turb()

    INTEGER :: i, j, k, ix, iy, iz
    REAL*8 :: rho


    DO k = 1-nb, nz+nb
        DO j = 1-nb, ny+nb
            DO i = 1-nb, nx+nb
                rho0(i,j,k) = q_3d(i,j,k,1)
            END DO
        END DO
    END DO
        
    !***********************************************************************************************************************************
    ! Compute kinetic energy equation quadratic co-efficients and average momentum contrained in the forcing. Will need those later to  
    ! normalize for constant turbulent kinetic energy input and make the force inject zero net momentum.
    !************************************************************************************************************************************
    quad_term1 = 0.d0
    quad_term2 = 0.d0
    momx_avg = 0.d0
    momy_avg = 0.d0
    momz_avg = 0.d0
        
    !compute quadratic co-efficients and average momentum injected by the forcing velocity 
    DO k = 1, nz
        DO j = 1, ny
            DO i = 1, nx
                
                ix = i - 1 - nx/2 
                iy = j - 1 - ny/2
                iz = k - 1 - nz/2

                rho = q_3d(i,j,k,1)

                ! volume integrated 1/2 rho del_v (dot) del_v
                quad_term1 = quad_term1 + rho * (delv_re(ix,iy,iz,1)**2 + delv_re(ix,iy,iz,2)**2 + delv_re(ix,iy,iz,3)**2) 
                    
                ! volume integrated del_v (dot) x-momentum
                quad_term2 = quad_term2 + delv_re(ix,iy,iz,1)*q_3d(i,j,k,2) + delv_re(ix,iy,iz,2)*q_3d(i,j,k,3) + delv_re(ix,iy,iz,3)*q_3d(i,j,k,4)

                ! average momentum
                momx_avg = momx_avg + rho * delv_re(ix,iy,iz,1) 
                momy_avg = momy_avg + rho * delv_re(ix,iy,iz,2) 
                momz_avg = momz_avg + rho * delv_re(ix,iy,iz,3) 

            END DO
        END DO
    END DO
        
    quad_term1 = 0.5d0 * quad_term1
    quad_term2 = quad_term2
       
 
    ! perform non-blocking MPI reduction to get sum of norm quadratic equation coefficients over all processes
    CALL start_allreduce(quad_term1, quad_term1_sum, quad_term1_req, 2)
    CALL start_allreduce(quad_term2, quad_term2_sum, quad_term2_req, 2)
    CALL start_allreduce(momx_avg, momx_avg_sum, momx_req, 2)
    CALL start_allreduce(momy_avg, momy_avg_sum, momy_req, 2)
    CALL start_allreduce(momz_avg, momz_avg_sum, momz_req, 2)
  

END SUBROUTINE get_rho_turb



! This subroutine generates a gaussian random force for driving the turbulence
SUBROUTINE generate_rand_force()

    REAL*8 :: delvk_re(-k_max:k_max, -k_max:k_max, -k_max:k_max, 3)
    REAL*8 :: delvk_im(-k_max:k_max, -k_max:k_max, -k_max:k_max, 3)
    REAL*8 :: dft_buffer_re(-k_max:k_max-1),dft_buffer_im(-k_max:k_max-1),dft_buffer_out_re(0:nnmax-1+2*nb),dft_buffer_out_im(0:nnmax-1+2*nb)   
    REAL*8 :: dvk_re(3), dvk_im(3), amp(3), phase(3)
    REAL*8 :: dvk_re_comp(3), dvk_re_sol(3), dvk_im_comp(3),dvk_im_sol(3) 
    REAL*8 :: rn(6), Pk, k_pk, kx0, ky0, kz0, k
    REAL*8 :: dotk_re, dotk_im, comp_frac, sol_frac, const_sol, const_comp, idft_pre
    INTEGER :: kx, ky, kz, i, ix, iy, iz
    REAL*8 :: t1, t2
   
    INTEGER, SAVE :: sgn = -1
    REAL*8 :: vx_amp, vy_amp, vx_ph, vy_ph, x, y, z
    INTEGER :: ixx, iyy, izz

    t1 = MPI_WTIME()

    kx0 = TWOPI / Lx
    ky0 = TWOPI / Ly
    kz0 = TWOPI / Lz 
    k_pk = TWOPI / MAX(Lx,Ly,Lz) 
    
    
    !*************************************************************
    ! Set fractional amounts for compressive and solenoidal forces
    !*************************************************************
  
    ! Pure solenoidal driving
    sol_frac = 1.d0
  
    ! Pure compressive driving
    !sol_frac = 0.d0
    
    ! Half-n-half
    !sol_frac = 0.5d0

    ! make sure 0 < sol_frac < 1
    !sol_frac = MAX(0.d0, MIN(1.d0, sol_frac))

    const_sol = SQRT(sol_frac)
    const_comp = SQRT(2.d0*(1.d0 - sol_frac)) 
    
    delv_rms = 0.d0
    
    !*******************************************
    ! Compute fourier components of random force
    !*******************************************
    DO kz = 0, k_max-1
        DO ky = 0, k_max-1
            DO kx = 0, k_max-1

                k = SQRT( (kx*kx0)**2 + (ky*ky0)**2  + (kz*kz0)**2 )
                   
                ! gaussian random field: power spectrum 
                Pk = (k**6) * EXP(-8.d0*k/k_pk)
                    
                ! generate six random numbers (one for amplitude and one for phase, for each velocity component)
                rn(1) = rand_taus()
                rn(2) = rand_taus()
                rn(3) = rand_taus()
                rn(4) = rand_taus()
                rn(5) = rand_taus()
                rn(6) = rand_taus()
                                  
                ! compute amplitudes and phases
                amp(1)   = SQRT(-2.d0 * Pk * LOG(rn(1)))
                amp(2)   = SQRT(-2.d0 * Pk * LOG(rn(2)))
                amp(3)   = SQRT(-2.d0 * Pk * LOG(rn(3)))
                phase(1) = TWOPI * rn(4) 
                phase(2) = TWOPI * rn(5) 
                phase(3) = TWOPI * rn(6) 
                     
                ! compute real and imaginary parts of velocity fourier component
                dvk_re(1) = amp(1) * COS(phase(1))
                dvk_re(2) = amp(2) * COS(phase(2))
                dvk_re(3) = amp(3) * COS(phase(3))
                dvk_im(1) = amp(1) * SIN(phase(1))
                dvk_im(2) = amp(2) * SIN(phase(2))
                dvk_im(3) = amp(3) * SIN(phase(3))
                   
                ! dot product with k_hat unit vector
                dotk_re = khat(kx,ky,kz,1) * dvk_re(1) + khat(kx,ky,kz,2) * dvk_re(2) + khat(kx,ky,kz,3) * dvk_re(3)
                dotk_im = khat(kx,ky,kz,1) * dvk_im(1) + khat(kx,ky,kz,2) * dvk_im(2) + khat(kx,ky,kz,3) * dvk_im(3)
                   
                ! decompose into compressive and solenoidal parts (by projection along and across k-vector)
                dvk_re_comp(1:3) = khat(kx,ky,kz,1:3) * dotk_re 
                dvk_im_comp(1:3) = khat(kx,ky,kz,1:3) * dotk_im
                                  
                dvk_re_sol(1:3) = dvk_re(1:3) - dvk_re_comp(1:3)               
                dvk_im_sol(1:3) = dvk_im(1:3) - dvk_im_comp(1:3)               
                               
                ! combine compressive and solenoidal components according to the prescribed fractions               
                delvk_re(kx,ky,kz,1:3) = const_sol * dvk_re_sol(1:3) + const_comp * dvk_re_comp(1:3)
                delvk_im(kx,ky,kz,1:3) = const_sol * dvk_im_sol(1:3) + const_comp * dvk_im_comp(1:3)
                   
                ! even (odd) symmetry in real (imaginary) part due to the driving force being real-valued function
                delvk_re(-kx,-ky,-kz,1:3) = delvk_re(kx,ky,kz,1:3)
                delvk_im(-kx,-ky,-kz,1:3) = -delvk_im(kx,ky,kz,1:3)
                   
                !delv_rms = delv_rms + 2.d0 * ( delvk_re(kx,ky,kz,1)**2 + delvk_im(kx,ky,kz,1)**2 + delvk_re(kx,ky,kz,2)**2 + &
                !           delvk_im(kx,ky,kz,2)**2 + delvk_re(kx,ky,kz,3)**2 + delvk_im(kx,ky,kz,3)**2 )                   
                   
            END DO
        END DO
    END DO

    !delv_rms = SQRT(delv_rms/DBLE(nx*ny*nz)) 

    !******************************************************************************************************************
    ! Do inverse fourier transform to obtain force in configuration space (via 1D DFTs along x,y and z directions)
    !******************************************************************************************************************

    delv_re = 0.d0
    delv_im = 0.d0
    idft_pre = 1.d0 / DBLE(nranks_x*nx * nranks_y*ny * nranks_z*nz) 
        
    ! DFT in x-direction
     DO i = 1, 3
        DO kz = -k_max, k_max-1
            DO ky = -k_max, k_max-1
                            
                ! copy strips-along x into 1d buffer    
                DO kx = -k_max, k_max-1
                    dft_buffer_re(kx) = delvk_re(kx,ky,kz,i)
                    dft_buffer_im(kx) = delvk_im(kx,ky,kz,i)
                END DO
    
                ! perform 1D inverse DFT 
                CALL direct_idft_1d(dft_buffer_re,dft_buffer_im,dft_buffer_out_re,dft_buffer_out_im,1-nb+my_coord(1)*nx,nx+nb+my_coord(1)*nx,nx*nranks_x)
    
                ! copy back into delv array 
                DO ix = -nb-nx/2, nb-1+nx/2
                    delv_re(ix,ky,kz,i) = dft_buffer_out_re(ix+nb+nx/2)  
                    delv_im(ix,ky,kz,i) = dft_buffer_out_im(ix+nb+nx/2)  
                END DO

            END DO       
        END DO
    END DO
    

    ! DFT in y-direction
    DO i = 1, 3
        DO kz = -k_max, k_max-1
            DO ix = -nb-nx/2, nb-1+nx/2
                            
                ! copy strips-along y into 1d buffer    
                DO ky = -k_max, k_max-1
                    dft_buffer_re(ky) = delv_re(ix,ky,kz,i)
                    dft_buffer_im(ky) = delv_im(ix,ky,kz,i)
                END DO
    
                ! perform 1D inverse DFT 
                CALL direct_idft_1d(dft_buffer_re,dft_buffer_im,dft_buffer_out_re,dft_buffer_out_im,1-nb+my_coord(2)*ny,ny+nb+my_coord(2)*ny,ny*nranks_y)
    
                ! copy back into delvk array 
                DO iy = -nb-ny/2, nb-1+ny/2
                    delv_re(ix,iy,kz,i) = dft_buffer_out_re(iy+nb+ny/2)  
                    delv_im(ix,iy,kz,i) = dft_buffer_out_im(iy+nb+ny/2)  
                END DO
                
            
            END DO
        END DO
    END DO
    
    ! DFT in z-direction
    DO i = 1, 3  
        DO iy = -nb-ny/2, nb-1+ny/2 
            DO ix = -nb-nx/2, nb-1+nx/2
                            
                ! copy strips-along z into 1d buffer    
                DO kz = -k_max, k_max-1
                    dft_buffer_re(kz) = delv_re(ix,iy,kz,i)
                    dft_buffer_im(kz) = delv_im(ix,iy,kz,i)
                END DO
    
                ! perform 1D inverse DFT 
                CALL direct_idft_1d(dft_buffer_re,dft_buffer_im,dft_buffer_out_re,dft_buffer_out_im,1-nb+my_coord(3)*nz,nz+nb+my_coord(3)*nz,nz*nranks_z)
    
                DO iz = -nb-nz/2, nb-1+nz/2
                    delv_re(ix,iy,iz,i) = dft_buffer_out_re(iz+nb+nz/2) * idft_pre  ! multiply by the IDFT prefactor of 1/nx*ny*nz
                END DO
                
            END DO        
        END DO
    END DO


    t2 = MPI_WTIME()
    t_force = t_force + t2-t1
    
    
END SUBROUTINE generate_rand_force


! This subroutine adds a random force to the fluid
SUBROUTINE add_random_force(dt, delv, delv_max)

    REAL*8, INTENT(IN) :: dt
    REAL*8, INTENT(IN) :: delv(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,3)
    REAL*8, INTENT(INOUT) :: delv_max
    REAL*8 :: rhoE, tmp1, tmp2, constE, delmomx, delmomy, delmomz
    INTEGER :: i, j, k
      

    ! wait till MPI allreduce completes
    CALL end_allreduce(quad_term1_req)
    CALL end_allreduce(quad_term2_req)
    CALL end_allreduce(momx_req)
    CALL end_allreduce(momy_req)
    CALL end_allreduce(momz_req)

    ! If only one MPI rank, then copy the original value into sum
    IF(nranks_x*nranks_y*nranks_z .EQ. 1) THEN
        quad_term1_sum = quad_term1
        quad_term2_sum = quad_term2
        momx_avg_sum = momx_avg
        momy_avg_sum = momy_avg
        momz_avg_sum = momz_avg
    END IF
            
    IF(quad_term1_sum .LE. 1.d-20) THEN
        PRINT*,'ERROR! Bad value in turbulence driver constant energy normalization!'
        PRINT*,'Myrank = ',myrank
        PRINT*,'quad_term1_sum, quad_term2_sum = ',quad_term1_sum, quad_term2_sum
        STOP
    END IF

    momx_avg_sum = momx_avg_sum * (my_dx**3)
    momy_avg_sum = momy_avg_sum * (my_dx**3)
    momz_avg_sum = momz_avg_sum * (my_dx**3)

    ! compute the constant-energy-input normalization constant (solve the kinetic energy quadratic equation and take the larger root)
    tmp1 = Edot * dt_turb / (my_dx**3)        
    tmp2 = SQRT(quad_term2_sum**2 + (4.d0 * tmp1 * quad_term1_sum))
    constE = 0.5d0 * (-quad_term2_sum + tmp2) / quad_term1_sum


    DO k = 1, nz
        DO j = 1, ny
            DO i = 1, nx

                !subtract off non-zero average momentum (so that no net momentum is injected) and normalize for constant energy input
                delmomx = constE * (rho0(i,j,k)*delv(i,j,k,1) - momx_avg_sum)
                delmomy = constE * (rho0(i,j,k)*delv(i,j,k,2) - momy_avg_sum)
                delmomz = constE * (rho0(i,j,k)*delv(i,j,k,3) - momz_avg_sum)
                
                q_3d(i,j,k,2) = q_3d(i,j,k,2) + delmomx                     
                q_3d(i,j,k,3) = q_3d(i,j,k,3) + delmomy                    
                q_3d(i,j,k,4) = q_3d(i,j,k,4) + delmomz                      

                delv_max = MAX(delv_max, ABS(delmomx)/ rho0(i,j,k), ABS(delmomy)/ rho0(i,j,k), ABS(delmomz)/ rho0(i,j,k)) 

            END DO
        END DO
    END DO
 

END SUBROUTINE add_random_force


! This subroutine performs a 1D inverse discrete fourier transform via direct summation
SUBROUTINE direct_idft_1d(in_re, in_im, out_re, out_im, ilow, ihi, nnx)

    REAL(8), INTENT(IN) :: in_re(-k_max:k_max-1), in_im(-k_max:k_max-1)    
    REAL(8), INTENT(OUT) :: out_re(0:nnmax-1+2*nb), out_im(0:nnmax-1+2*nb)   
    INTEGER, INTENT(IN) :: ilow, ihi, nnx
    REAL(8) :: sgn = -1.d0     ! DFT for sgn = 1 and Inverse DFT for sgn = -1
    INTEGER :: i, j
    REAL(8) :: theta, theta0
    
    theta0 =  TWOPI / DBLE(nnx)
    
    ! clear output arrays
    out_re = 0.d0
    out_im = 0.d0
        
        
    DO i = ilow, ihi
        DO j = -(k_max-1), k_max-1 !  0, nx-1
            theta = theta0 * i * j
            out_re(i-ilow) = out_re(i-ilow) + (in_re(j) * DCOS(theta) + sgn * in_im(j) * DSIN(theta)) 
            out_im(i-ilow) = out_im(i-ilow) + (in_im(j) * DCOS(theta) - sgn * in_re(j) * DSIN(theta))
        END DO
    END DO   
  
  
END SUBROUTINE direct_idft_1d



END MODULE turb_driver_mod