PROGRAM DFT_3d


IMPLICIT NONE


INTEGER, PARAMETER :: nx = 32
INTEGER, PARAMETER :: ny = nx
INTEGER, PARAMETER :: nz = nx
INTEGER, PARAMETER :: k_max = nx/2 
REAL(8), PARAMETER :: TWOPI = 8.d0*ATAN(1.d0)

REAL*8, ALLOCATABLE :: delv_re(:,:,:,:), delv_im(:,:,:,:), delvk_re(:,:,:,:), delvk_im(:,:,:,:)
REAL*8 :: khat(1:k_max-1,1:k_max-1,1:k_max-1,3)    
REAL*8 :: Lx, Ly, Lz, dx
REAL*8 :: K
INTEGER :: ix, iy, iz, kx, ky, kz


ALLOCATE(delvk_re(-nx/2:-1+nx/2,-ny/2:-1+ny/2,-nz/2:-1+nz/2,3),delvk_im(-nx/2:-1+nx/2,-ny/2:-1+ny/2,-nz/2:-1+nz/2,3))
ALLOCATE(delv_re(-nx/2:-1+nx/2,-ny/2:-1+ny/2,-nz/2:-1+nz/2,3),delv_im(-nx/2:-1+nx/2,-ny/2:-1+ny/2,-nz/2:-1+nz/2,3))

! grid size
Lx = 1.d0
Ly = 1.d0
Lz = 1.d0

! precomupte k-space unit vectors
kx0 = TWOPI / Lx
ky0 = TWOPI / Ly
kz0 = TWOPI / Lz 
    
khat = 0.d0
DO kz = 1, k_max-1
    DO ky = 1, k_max-1
        DO kx = 1, k_max-1 
            k = SQRT( (kx*kx0)**2 + (ky*ky0)**2  + (kz*kz0)**2 )
            khat(kx,ky,kz,1) = kx0 * kx / k
            khat(kx,ky,kz,2) = ky0 * ky / k
            khat(kx,ky,kz,3) = kz0 * kz / k
        END DO
    END DO
END DO


! set fourier components

delvk_re = 0.d0
delvk_im = 0.d0

DO kz = 1, k_max -1 
    DO ky = 1, k_max -1 
        DO kx = 1, k_max -1 

            IF(kx .EQ. 1 .AND. ky .EQ. 0 .AND. kz .EQ. 0 ) THEN
                delvk_re(kx, ky, kz) = 1.d0
                delvk_im(kx, ky, kz) = 0.d0
            ELSE
                delvk_re(kx, ky, kz) = 0.d0
                delvk_im(kx, ky, kz) = 0.d0
            END IF
            
            
            ! apply symmetry constraint from function being real in config. space
            delvk_re(-kx, -ky, -kz) = delvk_re(kx, ky, kz)
            delvk_im(-kx, -ky, -kz) = -delvk_im(kx, ky, kz)

        END DO
    END DO
END DO


! now compute the inverse DFT



CONTAINS


SUBROUTINE dft_3d()


    REAL*8 :: dft_buffer_re(-k_max:k_max-1),dft_buffer_im(-k_max:k_max-1),dft_buffer_out_re(1:nx),dft_buffer_out_im(1:nx)   
    INTEGER :: kx, ky, kz, i, ix, iy, iz
  
      
    !******************************************************************************************************************
    ! Do inverse fourier transform to obtain force in configuration space (separate 1D DFTs along x,y and z directions)
    !******************************************************************************************************************

    delv_re = 0.d0
    delv_im = 0.d0
    idft_pre = 1.d0 / DBLE(nx*ny *nz) 
    
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
                CALL direct_idft_1d(dft_buffer_re,dft_buffer_im,dft_buffer_out_re,dft_buffer_out_im)
    
                ! copy back into delv array 
                DO ix = -nx/2, -1+nx/2
                    delv_re(ix,ky,kz,i) = dft_buffer_out_re(ix+nx/2)  
                    delv_im(ix,ky,kz,i) = dft_buffer_out_im(ix+nx/2)  
                END DO

            END DO       
        END DO
    END DO
    

    ! DFT in y-direction
    DO i = 1, 3
        DO kz = -k_max, k_max-1
            DO ix = -nx/2, -1+nx/2
                            
                ! copy strips-along y into 1d buffer    
                DO ky = -k_max, k_max-1
                    dft_buffer_re(ky) = delv_re(ix,ky,kz,i)
                    dft_buffer_im(ky) = delv_im(ix,ky,kz,i)
                END DO
    
                ! perform 1D inverse DFT 
                CALL direct_idft_1d(dft_buffer_re,dft_buffer_im,dft_buffer_out_re,dft_buffer_out_im)
    
                ! copy back into delvk array 
                DO iy = -ny/2, -1+ny/2
                    delv_re(ix,iy,kz,i) = dft_buffer_out_re(iy+ny/2)  
                    delv_im(ix,iy,kz,i) = dft_buffer_out_im(iy+ny/2)  
                END DO
                
            
            END DO
        END DO
    END DO
    
    ! DFT in z-direction. (also apply the IDFT prefactor of 1/nx*ny*nz)
    DO i = 1, 3  
        DO iy = -ny/2, -1+ny/2 
            DO ix = -nx/2, -1+nx/2
                            
                ! copy strips-along z into 1d buffer    
                DO kz = -k_max, k_max-1
                    dft_buffer_re(kz) = delv_re(ix,iy,kz,i)
                    dft_buffer_im(kz) = delv_im(ix,iy,kz,i)
                END DO
    
                ! perform 1D inverse DFT 
                CALL direct_idft_1d(dft_buffer_re,dft_buffer_im,dft_buffer_out_re,dft_buffer_out_im)
    
                DO iz = -nz/2, -1+nz/2
                    delv_re(ix,iy,iz,i) = dft_buffer_out_re(iz+nz/2) * idft_pre  
                END DO
                
            END DO        
        END DO
    END DO

END SUBROUTINE dft_3d


! This subroutine performs a 1D inverse discrete fourier transform via direct summation
SUBROUTINE direct_idft_1d(in_re, in_im, out_re, out_im)

    REAL(8), INTENT(IN) :: in_re(-k_max:k_max-1), in_im(-k_max:k_max-1)    
    REAL(8), INTENT(OUT) :: out_re(1:nx), out_im(1:nx)   
    INTEGER :: i, j
    REAL(8) :: theta, theta0
    
    theta0 =  TWOPI / DBLE(nx)
    
    ! clear output arrays
    out_re = 0.d0
    out_im = 0.d0
        
    ilow = 1    
    ihi = nx
    
    DO i = ilow, ihi
        DO j = -kmax, k_max-1 
            theta = theta0 * i * j
            out_re(i-ilow) = out_re(i-ilow) + (in_re(j) * DCOS(theta) - in_im(j) * DSIN(theta)) 
            out_im(i-ilow) = out_im(i-ilow) + (in_im(j) * DCOS(theta) + in_re(j) * DSIN(theta))
        END DO
    END DO   
  
  

END SUBROUTINE direct_idft_1d





END PROGRAM DFT_3d