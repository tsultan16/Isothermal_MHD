MODULE grid_data_mod

USE constants_mod


IMPLICIT NONE	

 
REAL*8, ALLOCATABLE :: q_3d(:,:,:,:), &         ! vector of conserved variables
                       qpass_3d(:,:,:,:), &     ! vector of passive variables
                       bface_3d(:,:,:,:), &     ! cell face magnetic field                      
                       emf_corner(:,:,:,:),    & ! cell-corner CT electric field
                       file_buffer_dble(:), &    ! parallel io buffer (double precision)
                       delv_re(:,:,:,:), delv_im(:,:,:,:) ! turbulence driving force   

REAL*4, ALLOCATABLE :: file_buffer_sngl(:)  ! parallel io buffer (single precision)

REAL*8 :: my_dt, my_dx                            

REAL*8  :: ttot = 0.d0, tsolve = 0.d0, tfile = 0.d0, tcom = 0.d0, trma = 0.d0, tred = 0.d0, &
           tsolve_1, tsolve_2, tsolve_3, tsolve_4, tfile_1, tfile_2, t1, t2, tcom_1, tcom_2, tcom_3, tcom_4, &
           tcom_5, tcom_6, tcom_7, tcom_8, tcom_9, tcom_10, tfile_3, tfile_4, tfile_5, tfile_6   
           
INTEGER :: s1, s2, s3  ! RNG seed       
        
INTEGER :: tsteps, dump_count        
           
REAL*8 :: t_force = 0.d0           
           
CONTAINS


SUBROUTINE create_grid_arrays()

    INTEGER :: total_grid_memory


    ALLOCATE(q_3d(1-nb:nx+nb, 1-nb:ny+nb,1-nb:nz+nb,7))
    ALLOCATE(bface_3d(1-nb:nx+nb, 1-nb:ny+nb,1-nb:nz+nb,3))
    ALLOCATE(emf_corner(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,3))

    IF(npass .GT. 0) ALLOCATE(qpass_3d(1-nb+2:nx+nb-2, 1-nb+2:ny+nb-2,1-nb+2:nz+nb-2,npass))
    IF(drive_turbulence) ALLOCATE(delv_re(-nb-nx/2:nb-1+nx/2,-nb-ny/2:nb-1+ny/2,-nb-nz/2:nb-1+nz/2,3),delv_im(-nb-nx/2:nb-1+nx/2,-nb-ny/2:nb-1+ny/2,-nb-nz/2:nb-1+nz/2,3))
    IF(parallel_io .AND. double_prec) ALLOCATE(file_buffer_dble(parallel_filesize))
    IF(parallel_io .AND. (.NOT. double_prec)) ALLOCATE(file_buffer_sngl(parallel_filesize))
    
    q_3d = 0.d0
    IF(npass .GT. 0) qpass_3d = 0.d0
    bface_3d = 0.d0
    
    total_grid_memory = SIZEOF(q_3d)  + SIZEOF(bface_3d) + SIZEOF(emf_corner)
    IF(npass .GT. 0) total_grid_memory = total_grid_memory + SIZEOF(qpass_3d)
    IF(parallel_io .AND. double_prec) total_grid_memory = total_grid_memory + SIZEOF(file_buffer_dble)
    IF(parallel_io .AND. (.NOT. double_prec)) total_grid_memory = total_grid_memory + SIZEOF(file_buffer_sngl)
    
    IF(myrank .EQ. 0) THEN
    
    PRINT*,''
    PRINT*,'#####################################################################'
    PRINT*,'Total Memory for grid array storage per MPI rank (Mb) = ',total_grid_memory*1e-6                        
    PRINT*,'#####################################################################'
    PRINT*,''

    END IF
    
    
END SUBROUTINE create_grid_arrays


SUBROUTINE destroy_grid_arrays()

    IF(myrank .EQ. 0) PRINT*,'Deallocating grid arrays.'
    
    DEALLOCATE(q_3d, bface_3d, emf_corner)
    IF(npass .GT. 0) DEALLOCATE(qpass_3d) 
    IF(drive_turbulence) DEALLOCATE(delv_re, delv_im)
    IF(parallel_io .AND. double_prec) DEALLOCATE(file_buffer_dble)
    IF(parallel_io .AND. (.NOT. double_prec)) DEALLOCATE(file_buffer_sngl)

END SUBROUTINE destroy_grid_arrays 



END MODULE grid_data_mod