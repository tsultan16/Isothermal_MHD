MODULE grid_data_mod

USE constants_mod


IMPLICIT NONE	

 
REAL*8, ALLOCATABLE :: q(:,:)    ! vector of conserved variables
                       
                            
REAL*8 :: Lx         
REAL*8 :: xmin, xmax
INTEGER :: i_sim



! main MPI variables
!INTEGER :: comm2d, ierr, ndim, myrank, numprocs(1), dims(2), mycoord(2), req(4), &
!           coltype, source, tag, destination

!LOGICAL :: isperiodic(2), reorder 

! MPI buffers
!INTEGER :: neighbor_rank(8), field_buffer_size_x, field_buffer_size_y, particle_buffer_size
!REAL*8, ALLOCATABLE :: particle_buffer_in(:), &
!                       field_buffer_xm(:), field_buffer_xp(:), field_buffer_ym(:), field_buffer_yp(:)
                       
!INTEGER :: xlow, ylow, zlow, nb_xlow, nb_xhi, nb_ylow, nb_yhi 
!INTEGER :: bxlow, bylow, bzlow, bxhi, byhi, bzhi, &
!           exlow, eylow, ezlow, exhi, eyhi, ezhi  
 
 
 
CONTAINS



SUBROUTINE create_grid_arrays()

    INTEGER :: total_grid_memory


	ALLOCATE(q(1-nb:nx+nb,7))
    
    q = 0.d0

    total_grid_memory = SIZEOF(q) 
    
    
    PRINT*,''
    PRINT*,'#####################################################################'
    PRINT*,'Total Memory for grid array storage (Mb) = ',total_grid_memory*1e-6                        
    PRINT*,'#####################################################################'
    PRINT*,''

END SUBROUTINE create_grid_arrays



SUBROUTINE destroy_grid_arrays()

    PRINT*,'Deallocating grid arrays.'
    
    DEALLOCATE(q)

END SUBROUTINE destroy_grid_arrays 



 
END MODULE grid_data_mod