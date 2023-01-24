MODULE grid_data_mod

USE constants_mod


IMPLICIT NONE	

 
REAL*8, ALLOCATABLE :: q_3d(:,:,:,:), &         ! vector of conserved variables
                       bface_3d(:,:,:,:), &     ! cell face magnetic field
                       file_buffer(:)           ! parallel io buffer
REAL*8 :: my_dt, my_dx                            
REAL*8 :: Lx         
REAL*8 :: xmin, xmax
INTEGER :: i_sim


 
CONTAINS



SUBROUTINE create_grid_arrays()

    INTEGER :: total_grid_memory


	ALLOCATE(q_3d(1-nb:nx+nb, 1-nb:ny+nb,1-nb:nz+nb,7))
	ALLOCATE(bface_3d(1-nb:nx+nb, 1-nb:ny+nb,1-nb:nz+nb,3))
    IF(parallel_io) ALLOCATE(file_buffer(parallel_filesize))
    
    q_3d = 0.d0
    bface_3d = 0.d0
    
    total_grid_memory = SIZEOF(q_3d)  + SIZEOF(bface_3d)
    IF(parallel_io) total_grid_memory = total_grid_memory + SIZEOF(file_buffer)
    
    IF(myrank .EQ. 0) THEN
    
    PRINT*,''
    PRINT*,'#####################################################################'
    PRINT*,'Total Memory for grid array storage (Mb) = ',total_grid_memory*1e-6                        
    PRINT*,'#####################################################################'
    PRINT*,''

    END IF
    
    
END SUBROUTINE create_grid_arrays



SUBROUTINE destroy_grid_arrays()

    IF(myrank .EQ. 0) PRINT*,'Deallocating grid arrays.'
    
    DEALLOCATE(q_3d, bface_3d)
    IF(parallel_io) DEALLOCATE(file_buffer)

END SUBROUTINE destroy_grid_arrays 



 
END MODULE grid_data_mod