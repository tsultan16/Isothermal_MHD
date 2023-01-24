MODULE mpi_domain_mod

USE constants_mod
USE grid_data_mod
USE ISO_C_BINDING  ! needed for access to C pointer
!USE MPI

IMPLICIT NONE

INCLUDE 'mpif.h'  ! MPI runtime library access (if using an older version of MPI, need to uncomment this line and get rid of "USE MPI")


! MPI Parameters and Local Variables
INTEGER, PARAMETER :: ndim  = 3 ! domain decomposition dimensions

INTEGER, PARAMETER :: maxvars = 7

INTEGER :: comm3d, dims(3), numprocs(1), req(6), status(MPI_STATUS_SIZE), np
LOGICAL :: isperiodic(3), reorder, neighbor_locked(26) 
INTEGER :: neighbor, ierr, neighbor_rank(26), neighbor_repeat_list(26)

!RMA Vars
INTEGER :: win, sizedouble, winsize
INTEGER(KIND=MPI_ADDRESS_KIND) :: disp(26)! 6 faces + 12 edges + 8 corners
INTEGER(KIND=MPI_ADDRESS_KIND) :: shuffled_disp(26)
INTEGER :: send_buffer_size(26)
INTEGER :: window
TYPE(C_PTR) :: rma_cmem  ! C pointer to RMA window memory block
REAL(8), POINTER :: A(:) ! Fortran pointer that will be associated with the RMA window C pointer

! MPI buffer
REAL(8), ALLOCATABLE :: buffer_out(:)

                        
INTEGER :: max_buffer_size


CONTAINS


SUBROUTINE initialize_MPI()

    IF(ndim .NE. 3) THEN
        PRINT*,'ERROR!! Need to set ndim = 3. Terminating program...'
        STOP
    END IF

    CALL MPI_INIT(ierr)

    ! create a 3d cartesian communicator
    CALL mpi_comm3d_init()

    CALL MPI_BARRIER(comm3d, ierr)   

    ! find out MPI ranks of all my neighbors
    CALL mpi_get_neighbor_ranks()

    ! initialize my RMA windows
    CALL rma_init()

    CALL MPI_BARRIER(comm3d, ierr)   

    IF(myrank .EQ. 0) PRINT*,'MPI initilization complete.'
     

END SUBROUTINE initialize_MPI


SUBROUTINE shut_off_MPI()

    CALL close_rma_window()

    CALL MPI_BARRIER(comm3d, ierr)

    CALL MPI_FINALIZE(ierr)

    DEALLOCATE(buffer_out)

END SUBROUTINE shut_off_MPI


SUBROUTINE mpi_comm3d_init()

    INTEGER :: ierr
    
    ! find out number of processors
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs(1), ierr)
    
    np = numprocs(1)
    
     IF(np .NE. nranks_x * nranks_y * nranks_z) THEN
        PRINT*,'ERROR! Need -np equal to ', nranks_x * nranks_y * nranks_z
        STOP
     END IF    
    
    ! create a cartesian communicator for our 3d domain decomposition
    dims(1) = nranks_x
    dims(2) = nranks_y
    dims(3) = nranks_z
    
    isperiodic(1:3) = .TRUE. ! periodic boundaries by default (Note: Not using periodic domain can make MPI communications siginficantly slower!!)

    reorder = .TRUE.
    
    CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndim, dims, isperiodic, reorder, comm3d, ierr)
    
    ! get rank of this processor
    CALL MPI_COMM_RANK(comm3d, myrank, ierr)

    IF(myrank .EQ. 0) PRINT*,'# of processors =', np

    ! get rank coordinates
    CALL MPI_CART_COORDS(comm3d, myrank, ndim, my_coord, ierr)

    !PRINT*,'my_coord =',my_coord

END SUBROUTINE mpi_comm3d_init


SUBROUTINE mpi_get_neighbor_ranks()

    INTEGER :: i, j, neighbor_coord(3), rank_dest     
    
    ! get neighbor coordinates
        
    ! x- face
    neighbor_coord(:) = (/ my_coord(1)-1, my_coord(2), my_coord(3) /)


    IF( ((neighbor_coord(1) .GT. nranks_x-1) .OR.(neighbor_coord(1) .LT. 0) ) .AND. (.NOT. isperiodic(1))) THEN
        neighbor_rank(1) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(1) = rank_dest
    END IF
     
       
    ! x+ face
    neighbor_coord(:) = (/ my_coord(1)+1, my_coord(2), my_coord(3) /)

    IF( ((neighbor_coord(1) .GT. nranks_x-1) .OR.(neighbor_coord(1) .LT. 0) ) .AND. (.NOT. isperiodic(1))) THEN
        neighbor_rank(2) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(2) = rank_dest
    END IF   
    
    ! y- face
    neighbor_coord(:) = (/ my_coord(1), my_coord(2)-1, my_coord(3) /)

    IF( ((neighbor_coord(2) .GT. nranks_y-1) .OR.(neighbor_coord(2) .LT. 0) ) .AND. (.NOT. isperiodic(2))) THEN
        neighbor_rank(3) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(3) = rank_dest
    END IF
       
    ! y+ face
    neighbor_coord(:) = (/ my_coord(1), my_coord(2)+1, my_coord(3) /)


    IF( ((neighbor_coord(2) .GT. nranks_y-1) .OR.(neighbor_coord(2) .LT. 0) ) .AND. (.NOT. isperiodic(2))) THEN
        neighbor_rank(4) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(4) = rank_dest
    END IF
    
    ! z- face
    neighbor_coord(:) = (/ my_coord(1), my_coord(2), my_coord(3)-1 /)


    IF( ((neighbor_coord(3) .GT. nranks_z-1) .OR.(neighbor_coord(3) .LT. 0) ) .AND. (.NOT. isperiodic(3))) THEN
        neighbor_rank(5) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(5) = rank_dest
    END IF
    
    ! z+ face
    neighbor_coord(:) = (/ my_coord(1), my_coord(2), my_coord(3)+1 /)


    IF( ((neighbor_coord(3) .GT. nranks_z-1) .OR.(neighbor_coord(3) .LT. 0) ) .AND. (.NOT. isperiodic(3))) THEN
        neighbor_rank(6) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(6) = rank_dest
    END IF   
    
    ! x-edge 1
    neighbor_coord(:) = (/ my_coord(1), my_coord(2)-1, my_coord(3)-1 /)

    IF( ((neighbor_coord(2) .LT. 0) .AND. (.NOT. isperiodic(2))) .OR. &
        ((neighbor_coord(3) .LT. 0) .AND. (.NOT. isperiodic(3)))  ) THEN
        
        neighbor_rank(7) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(7) = rank_dest
    END IF
    
    
    ! x-edge 2
    neighbor_coord(:) = (/ my_coord(1), my_coord(2)+1, my_coord(3)-1 /)

    IF( ((neighbor_coord(2) .GT. nranks_y-1) .AND. (.NOT. isperiodic(2))) .OR. &
        ((neighbor_coord(3) .LT. 0) .AND. (.NOT. isperiodic(3)))  ) THEN
        neighbor_rank(8) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(8) = rank_dest
    END IF
    
    ! x-edge 3
    neighbor_coord(:) = (/ my_coord(1), my_coord(2)-1, my_coord(3)+1 /)

    IF( ((neighbor_coord(2) .LT. 0) .AND. (.NOT. isperiodic(2))) .OR. &
        ((neighbor_coord(3) .GT. nranks_z-1) .AND. (.NOT. isperiodic(3)))  ) THEN
        neighbor_rank(9) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(9) = rank_dest
    END IF
    
    ! x-edge 4
    neighbor_coord(:) = (/ my_coord(1), my_coord(2)+1, my_coord(3)+1 /)
    
    IF( ((neighbor_coord(2) .GT. nranks_y-1) .AND. (.NOT. isperiodic(2))) .OR. &
        ((neighbor_coord(3) .GT. nranks_z-1) .AND. (.NOT. isperiodic(3)))  ) THEN
        neighbor_rank(10) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(10) = rank_dest
    END IF
    
        
     ! y-edge 1
    neighbor_coord(:) = (/ my_coord(1)-1, my_coord(2), my_coord(3)-1 /)


    IF( ((neighbor_coord(1) .LT. 0) .AND. (.NOT. isperiodic(1))) .OR. &
        ((neighbor_coord(3) .LT. 0) .AND. (.NOT. isperiodic(3)))  ) THEN        


        neighbor_rank(11) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(11) = rank_dest
    END IF
    
    
    ! y-edge 2
    neighbor_coord(:) = (/ my_coord(1)+1, my_coord(2), my_coord(3)-1 /)

    IF( ((neighbor_coord(1) .GT. nranks_x-1) .AND. (.NOT. isperiodic(1))) .OR. &
        ((neighbor_coord(3) .LT. 0) .AND. (.NOT. isperiodic(3)))  ) THEN
        neighbor_rank(12) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(12) = rank_dest
    END IF
    
    ! y-edge 3
    neighbor_coord(:) = (/ my_coord(1)-1, my_coord(2), my_coord(3)+1 /)

    IF( ((neighbor_coord(1) .LT. 0) .AND. (.NOT. isperiodic(1))) .OR. &
        ((neighbor_coord(3) .GT. nranks_z-1) .AND. (.NOT. isperiodic(3)))  ) THEN
        neighbor_rank(13) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(13) = rank_dest
    END IF
    
    ! y-edge 4
    neighbor_coord(:) = (/ my_coord(1)+1, my_coord(2), my_coord(3)+1 /)

    IF( ((neighbor_coord(1) .GT. nranks_x-1) .AND. (.NOT. isperiodic(1))) .OR. &
        ((neighbor_coord(3) .GT. nranks_z-1) .AND. (.NOT. isperiodic(3)))  ) THEN
        neighbor_rank(14) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(14) = rank_dest
    END IF
    
     ! z-edge 1
    neighbor_coord(:) = (/ my_coord(1)-1, my_coord(2)-1, my_coord(3) /)

    IF( ((neighbor_coord(1) .LT. 0) .AND. (.NOT. isperiodic(1))) .OR. &
        ((neighbor_coord(2) .LT. 0) .AND. (.NOT. isperiodic(2)))  ) THEN
        neighbor_rank(15) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(15) = rank_dest
    END IF
    
    
    ! z-edge 2
    neighbor_coord(:) = (/ my_coord(1)+1, my_coord(2)-1, my_coord(3) /)

    IF( ((neighbor_coord(1) .GT. nranks_x-1) .AND. (.NOT. isperiodic(1))) .OR. &
        ((neighbor_coord(2) .LT. 0) .AND. (.NOT. isperiodic(2)))  ) THEN
        neighbor_rank(16) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(16) = rank_dest
    END IF
    
    ! z-edge 3
    neighbor_coord(:) = (/ my_coord(1)-1, my_coord(2)+1, my_coord(3) /)

    IF( ((neighbor_coord(1) .LT. 0) .AND. (.NOT. isperiodic(1))) .OR. &
        ((neighbor_coord(2) .GT. nranks_y-1) .AND. (.NOT. isperiodic(2)))  ) THEN
        neighbor_rank(17) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(17) = rank_dest
    END IF
    
    ! z-edge 4
    neighbor_coord(:) = (/ my_coord(1)+1, my_coord(2)+1, my_coord(3) /)

    IF( ((neighbor_coord(1) .GT. nranks_x-1) .AND. (.NOT. isperiodic(1))) .OR. &
        ((neighbor_coord(2) .GT. nranks_y-1) .AND. (.NOT. isperiodic(2)))  ) THEN
        neighbor_rank(18) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(18) = rank_dest
    END IF
    

    ! bottom SW corner
    neighbor_coord(:) = (/ my_coord(1)-1, my_coord(2)-1, my_coord(3)-1 /)

    IF( ((neighbor_coord(1) .LT. 0) .AND. (.NOT. isperiodic(1))) .OR. &
        ((neighbor_coord(2) .LT. 0) .AND. (.NOT. isperiodic(2))) .OR. &
        ((neighbor_coord(3) .LT. 0) .AND. (.NOT. isperiodic(3))) ) THEN
        neighbor_rank(19) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(19) = rank_dest
    END IF
    
    ! bottom SE corner
    neighbor_coord(:) = (/ my_coord(1)+1, my_coord(2)-1, my_coord(3)-1 /)

    IF( ((neighbor_coord(1) .GT. nranks_x-1) .AND. (.NOT. isperiodic(1))) .OR. &
        ((neighbor_coord(2) .LT. 0) .AND. (.NOT. isperiodic(2))) .OR. &
        ((neighbor_coord(3) .LT. 0) .AND. (.NOT. isperiodic(3))) ) THEN
        neighbor_rank(20) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(20) = rank_dest
    END IF
    
    
    ! bottom NW corner
    neighbor_coord(:) = (/ my_coord(1)-1, my_coord(2)+1, my_coord(3)-1 /)

    IF( ((neighbor_coord(1) .LT. 0) .AND. (.NOT. isperiodic(1))) .OR. &
        ((neighbor_coord(2) .GT. nranks_y-1) .AND. (.NOT. isperiodic(2))) .OR. &
        ((neighbor_coord(3) .LT. 0) .AND. (.NOT. isperiodic(3))) ) THEN
        neighbor_rank(21) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(21) = rank_dest
    END IF
    
    
    ! bottom NE corner
        neighbor_coord(:) = (/ my_coord(1)+1, my_coord(2)+1, my_coord(3)-1 /)

    IF( ((neighbor_coord(1) .GT. nranks_x-1) .AND. (.NOT. isperiodic(1))) .OR. &
        ((neighbor_coord(2) .GT. nranks_y-1) .AND. (.NOT. isperiodic(2))) .OR. &
        ((neighbor_coord(3) .LT. 0) .AND. (.NOT. isperiodic(3))) ) THEN
        neighbor_rank(22) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(22) = rank_dest
    END IF
    
    ! top SW corner
    neighbor_coord(:) = (/ my_coord(1)-1, my_coord(2)-1, my_coord(3)+1 /)

    IF( ((neighbor_coord(1) .LT. 0) .AND. (.NOT. isperiodic(1))) .OR. &
        ((neighbor_coord(2) .LT. 0) .AND. (.NOT. isperiodic(2))) .OR. &
        ((neighbor_coord(3) .GT. nranks_z-1) .AND. (.NOT. isperiodic(3))) ) THEN
        neighbor_rank(23) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(23) = rank_dest
    END IF
    
    ! top SE corner
    neighbor_coord(:) = (/ my_coord(1)+1, my_coord(2)-1, my_coord(3)+1 /)

    IF( ((neighbor_coord(1) .GT. nranks_x-1) .AND. (.NOT. isperiodic(1))) .OR. &
        ((neighbor_coord(2) .LT. 0) .AND. (.NOT. isperiodic(2))) .OR. &
        ((neighbor_coord(3) .GT. nranks_z-1) .AND. (.NOT. isperiodic(3))) ) THEN
        neighbor_rank(24) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(24) = rank_dest
    END IF
    
    
    ! top NW corner
    neighbor_coord(:) = (/ my_coord(1)-1, my_coord(2)+1, my_coord(3)+1 /)

    IF( ((neighbor_coord(1) .LT. 0) .AND. (.NOT. isperiodic(1))) .OR. &
        ((neighbor_coord(2) .GT. nranks_y-1) .AND. (.NOT. isperiodic(2))) .OR. &
        ((neighbor_coord(3) .GT. nranks_z-1) .AND. (.NOT. isperiodic(3))) ) THEN
        neighbor_rank(25) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(25) = rank_dest
    END IF
    
    
    ! top NE corner
    neighbor_coord(:) = (/ my_coord(1)+1, my_coord(2)+1, my_coord(3)+1 /)

    IF( ((neighbor_coord(1) .GT. nranks_x-1) .AND. (.NOT. isperiodic(1))) .OR. &
        ((neighbor_coord(2) .GT. nranks_y-1) .AND. (.NOT. isperiodic(2))) .OR. &
        ((neighbor_coord(3) .GT. nranks_z-1) .AND. (.NOT. isperiodic(3))) ) THEN
        neighbor_rank(26) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(26) = rank_dest
    END IF
    
    
    !**************************************************************************************
    ! Check whether multiple neighbors share the same rank. Need to know this to make sure
    ! that a WIN_LOCK is not opened more than once.
    !**************************************************************************************
    DO i = 1, 26
        neighbor_repeat_list(i) = i
    END DO
    
    DO i = 1, 26
        DO j = 1, i-1 
        
            ! record array index of repeated neighbor rank in the neighbor_repeat_list array            
            IF(neighbor_rank(i) .EQ. neighbor_rank(j)) THEN  
                neighbor_repeat_list(i) = j
                EXIT
            END IF       
            
        END DO   
    END DO
    
    
   ! DO i = 0, np-1
    
   !     CALL MPI_BARRIER(comm3d, ierr)
    
   !     IF(myrank .EQ. i) THEN
        
   !         PRINT*,''
   !         PRINT*,'#######################################################################################'
   !         WRITE(*,FMT='(" Myrank = ",i2,", Neighbor Faces:  x- = ",i2,", x+ = ",i2,", y- = ",i2,", y+ = ",i2,", z- = ",i2,", z+ = ",i2)') &
   !         myrank, neighbor_rank(1) ,neighbor_rank(2) ,neighbor_rank(3) ,neighbor_rank(4) ,neighbor_rank(5) ,neighbor_rank(6) 
   !         WRITE(*,FMT='(" X-edge:  edge1 = ",i2,", edge2 = ",i2,", edge3 = ",i2,", edge4 = ",i2)') &
   !         neighbor_rank(7) ,neighbor_rank(8) ,neighbor_rank(9) ,neighbor_rank(10) 
   !         WRITE(*,FMT='(" Y-edge:  edge1 = ",i2,", edge2 = ",i2,", edge3 = ",i2,", edge4 = ",i2)') &
   !         neighbor_rank(11) ,neighbor_rank(12) ,neighbor_rank(13) ,neighbor_rank(14) 
   !         WRITE(*,FMT='(" Z-edge:  edge1 = ",i2,", edge2 = ",i2,", edge3 = ",i2,", edge4 = ",i2)') &
   !         neighbor_rank(15) ,neighbor_rank(16) ,neighbor_rank(17) ,neighbor_rank(18)
   !         WRITE(*,FMT='(" Bottom-corners:  SW = ",i2,", SE = ",i2,", NW = ",i2,", NE = ",i2)') &
   !         neighbor_rank(19) ,neighbor_rank(20) ,neighbor_rank(21) ,neighbor_rank(22)
   !         WRITE(*,FMT='(" Top-corners:  SW = ",i2,", SE = ",i2,", NW = ",i2,", NE = ",i2)') &
   !         neighbor_rank(23) ,neighbor_rank(24) ,neighbor_rank(25) ,neighbor_rank(26)
            
   !         PRINT*,'My coordinates =',my_coord
   !         PRINT*,'######################################################################################'
   !         PRINT*,''
        
   !     END IF

   !     CALL MPI_BARRIER(comm3d, ierr)
    
    !END DO
    
    
    
    
END SUBROUTINE mpi_get_neighbor_ranks



! Note: The RMA window is broken up into 26 chunks corresponding to ghost cells at the 
! * 6 cell face boundaries: x-, x+, y-, y+, z-, z+
! * 12 cell edge  boundaries, i.e. 4 edges along to each coordinate direction: x(1,2,3,4), y(1,2,3,4), z(1,2,3,4) 
! * 8 cell corners: bottom (SW, SE, NW, NE), top(SW, SE, NW, NE)
! The first two slots in each chunk are reserved for the data recv/sent signals from the neighbor corresponding to that boundary.
SUBROUTINE rma_init()

    INTEGER :: i, ierr, num, nvars, nb_pass
    REAL*8 :: double_var, total_mpi_memory
    
    nb_pass = 3 ! number of boundary cells for passive variables
    
    ! sizes of PUT buffers corresponding to each of the 26 boundaries 
    send_buffer_size(1:2)   = ny * nz * (nb * 10 + nb_pass * npass)
    send_buffer_size(3:4)   = nx * nz * (nb * 10 + nb_pass * npass) 
    send_buffer_size(5:6)   = ny * nx * (nb * 10 + nb_pass * npass)
    send_buffer_size(7:10)  = nx * nb * (nb * 10 + nb_pass * npass) 
    send_buffer_size(11:14) = ny * nb * (nb * 10 + nb_pass * npass)
    send_buffer_size(15:18) = nz * nb * (nb * 10 + nb_pass * npass)
    send_buffer_size(19:26) = nb * nb * (nb * 10 + nb_pass * npass)

    
    ! determine array element size in bytes
    CALL MPI_SIZEOF(double_var, sizedouble, ierr)
    !sizedouble = 8
    num = 2 * 26 + SUM(send_buffer_size)
    winsize =  num * sizedouble
    
    IF(myrank .EQ. 0) PRINT*,'RMA window size per rank (Mbytes) = ',winsize*1e-6 
    
    ! create local rma memory window and allocate memory block
    CALL MPI_WIN_ALLOCATE(winsize, sizedouble, MPI_INFO_NULL, comm3d, rma_cmem, win, ierr)
    CALL C_F_POINTER(rma_cmem, A, (/ num /))
    
    ! Allocate memory for out-buffer
    ALLOCATE(buffer_out(num))
    
    A(:) = 0.d0

    !***********************************************************
    ! set rma window displacement units for the different chunks
    ! We will use displacement unit arrays for q_3d and bface_3d
    !***********************************************************

    disp(1)  = 0
    DO i = 2, 26
        disp(i) = disp(i-1) + send_buffer_size(i-1) + 2
    END DO
    
    ! shuffled displacement units (for convenience later when doing PUTs into neighbor boundaries)  
    shuffled_disp = (/ disp(2),  disp(1),  disp(4),  disp(3), disp(6), disp(5),  &
                       disp(10), disp(9),  disp(8),  disp(7),   &
                       disp(14), disp(13), disp(12), disp(11),  & 
                       disp(18), disp(17), disp(16), disp(15),  &
                       disp(26), disp(25), disp(24), disp(23),  &
                       disp(22), disp(21), disp(20), disp(19)     /)



    total_mpi_memory = winsize + SIZEOF(buffer_out)

    IF(myrank .EQ. 0) THEN  
        PRINT*,''
        PRINT*,'#########################################################################'
        PRINT*,'Total Memory for MPI buffer array storage  per MPI rank (Mb) = ',total_mpi_memory*1e-6                        
        PRINT*,'#########################################################################'
        PRINT*,''
    END IF


    IF(myrank .EQ. 0) PRINT*,'RMA window initialization completed.'


END SUBROUTINE rma_init



! close the RMA window
SUBROUTINE close_rma_window()

    INTEGER :: ierr, i, j

    ! free rma window
    CALL MPI_WIN_FREE(win, ierr)
    
    NULLIFY(A)


END SUBROUTINE close_rma_window


! This subroutine posts a non-blocking MPI all_reduce to compute a global minimum quantity
! across all mpi sub-domains
SUBROUTINE start_allreduce(my_value, reduced_value, reduce_req, op)

    REAL(8), INTENT(IN)    :: my_value
    REAL(8), INTENT(INOUT) :: reduced_value
    INTEGER, INTENT(INOUT) :: reduce_req
    INTEGER, INTENT(IN) :: op
    INTEGER :: ierr

    
    ! post non-blocking MPI reduction
    IF(op .EQ. 1) THEN
        CALL MPI_IALLREDUCE(my_value, reduced_value, 1, MPI_DOUBLE_PRECISION, MPI_MIN, comm3d, reduce_req, ierr)
    ELSE IF(op .EQ. 2) THEN
        CALL MPI_IALLREDUCE(my_value, reduced_value, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, reduce_req, ierr)
    END IF

END SUBROUTINE start_allreduce


! This subroutine causes program to waits for the MPI all_reduce to complete 
SUBROUTINE end_allreduce(reduce_req)

    INTEGER, INTENT(INOUT) :: reduce_req
    INTEGER :: status(MPI_STATUS_SIZE), ierr
    
    ! wait until MPI reduction completes
    CALL MPI_WAIT(reduce_req, status, ierr)
   
   
END SUBROUTINE end_allreduce


! Top-level for initiating halo exchange routine
SUBROUTINE start_halo_exchange()

    INTEGER :: i, ierr

    ! Start RMA epoch by opening up WIN_fence so neighbors can access my RMA window 
    CALL MPI_WIN_FENCE(0, win, ierr)
        
    DO i = 1, 26
        IF((neighbor_rank(i) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(i) .NE. myrank)) THEN
             CALL send_bndry(i)
        END IF
    END DO
    

END SUBROUTINE start_halo_exchange



! Top-level for completing halo exchange routine
SUBROUTINE end_halo_exchange()

    INTEGER :: ierr

    
    ! close WIN_fence to end the RMA epoch
    CALL MPI_WIN_FENCE(0, win, ierr)

    ! now we can unback boundary data
    CALL unpack_bufferin()    
    

END SUBROUTINE end_halo_exchange



! put boundary data into neighbor's RMA window
SUBROUTINE send_bndry(bndry)

    INTEGER, INTENT(IN) :: bndry
    INTEGER :: ierr, buff_size
    INTEGER(KIND=MPI_ADDRESS_KIND) :: msg_disp

    msg_disp = 2 + shuffled_disp(bndry)    
    buff_size = send_buffer_size(bndry)  
  
    ! pack up boundary data into out-buffer 
    CALL packup_outbuffer(bndry, buffer_out, INT(msg_disp))
    
    ! PUT into neighbor's RMA window    
    CALL MPI_PUT(buffer_out(msg_disp), buff_size, MPI_DOUBLE_PRECISION, neighbor_rank(bndry), &
                 msg_disp, buff_size, MPI_DOUBLE_PRECISION, win, ierr)
   
    
END SUBROUTINE send_bndry
 

SUBROUTINE packup_outbuffer(bndry, buff, msg_disp)

    INTEGER, INTENT(IN) :: bndry, msg_disp
    REAL*8, INTENT(OUT) :: buff(:)
    INTEGER :: i, j, k, p, ix, nb_pass

    nb_pass = 3

    ix = msg_disp

     SELECT CASE(bndry)

        !********************
        ! Faces
        !********************

        CASE(1) ! x- boundary

        ! first load up the cell-center mhd state variables
        DO i = 1, nb          
            DO k = 1, nz
                DO j = 1, ny
                    DO p = 1, 7 
                    
                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix + 1
                    
                    END DO
                END DO
            END DO
        END DO
        
       ! now load up the cell face magnetic field (extra face on the bottom ends of the boundary)
       DO i = 1, nb          
            DO k = 1, nz
                DO j = 1, ny
                    DO p = 1, 3 
                    
                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix + 1
                    
                    END DO
                END DO
            END DO
       END DO
       
       IF(npass .GT. 0) THEN       
           DO i = 1, nb_pass          
                DO k = 1, nz
                    DO j = 1, ny
                        DO p = 1, npass 
                        
                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix + 1
                        
                        END DO
                    END DO
                END DO
            END DO        
        END IF  
          
        CASE(2) ! x+ boundary

        DO i = nx, nx-nb+1, -1
            DO k = 1, nz
                DO j = 1, ny
                    DO p = 1, 7 

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO
        
        DO i = nx, nx-nb+1, -1
            DO k = 1, nz
                DO j = 1, ny
                    DO p = 1, 3 

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO        
        
        IF(npass .GT. 0) THEN 
            DO i = nx, nx-nb_pass+1, -1
                DO k = 1, nz
                    DO j = 1, ny
                        DO p = 1, npass 

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO 
                    END DO
                END DO
            END DO                
        END IF
        
        CASE(3) ! y- boundary

        DO j = 1, nb
            DO k = 1, nz
                DO i = 1, nx
                    DO p = 1, 7 

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO

        DO j = 1, nb
            DO k = 1, nz
                DO i = 1, nx
                    DO p = 1, 3 

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO
        
        IF(npass .GT. 0) THEN
            DO j = 1, nb_pass
                DO k = 1, nz
                    DO i = 1, nx
                        DO p = 1, npass 

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO 
                    END DO
                END DO
            END DO
        END IF
        
        CASE(4) ! y+ boundary

        DO j = ny, ny-nb+1, -1
            DO k = 1, nz
                DO i = 1, nx
                    DO p = 1, 7

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO   
                END DO
            END DO
        END DO    

        DO j = ny, ny-nb+1, -1
            DO k = 1, nz
                DO i = 1, nx
                    DO p = 1, 3

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO   
                END DO
            END DO
        END DO    
        
        IF(npass .GT. 0) THEN
            DO j = ny, ny-nb_pass+1, -1
                DO k = 1, nz
                    DO i = 1, nx
                        DO p = 1, npass

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO   
                    END DO
                END DO
            END DO         
        END IF
      
        CASE(5) ! z- boundary

        DO k = 1, nb
            DO j = 1, ny
                DO i = 1, nx
                    DO p = 1, 7 

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO    

        DO k = 1, nb
            DO j = 1, ny
                DO i = 1, nx
                    DO p = 1, 3 

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO    
        
        IF(npass .GT. 0) THEN
            DO k = 1, nb_pass
                DO j = 1, ny
                    DO i = 1, nx
                        DO p = 1, npass 

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO 
                    END DO
                END DO
            END DO
        END IF
        
        CASE(6) ! z+ boundary

        DO k = nz, nz-nb+1, -1
            DO j = 1, ny
                DO i = 1, nx
                    DO p = 1, 7 

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO                
                END DO
            END DO
        END DO
    
        DO k = nz, nz-nb+1, -1
            DO j = 1, ny
                DO i = 1, nx
                    DO p = 1, 3 

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO                
                END DO
            END DO
        END DO

        IF(npass .GT. 0) THEN
            DO k = nz, nz-nb_pass+1, -1
                DO j = 1, ny
                    DO i = 1, nx
                        DO p = 1, npass 

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO                
                    END DO
                END DO
            END DO 
        END IF

        !********************
        ! Edges
        !********************

        CASE(7) ! x-edge 1

        DO k = 1, nb  
            DO j = 1, nb
                DO i = 1, nx
                    DO p = 1, 7 

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO    
        
        DO k = 1, nb  
            DO j = 1, nb
                DO i = 1, nx
                    DO p = 1, 3 

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO    
        
        IF(npass .GT. 0) THEN
            DO k = 1, nb_pass  
                DO j = 1, nb_pass
                    DO i = 1, nx
                        DO p = 1, npass 

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO 
                    END DO
                END DO
            END DO  
        END IF
        
        CASE(8) ! x-edge 2

        DO k = 1,nb  
            DO j = ny, ny-nb+1,-1
                DO i = 1, nx
                    DO p = 1, 7 

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO    

        DO k = 1,nb  
            DO j = ny, ny-nb+1,-1
                DO i = 1, nx
                    DO p = 1, 3 

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO    
        
        IF(npass .GT. 0) THEN
            DO k = 1,nb_pass  
                DO j = ny, ny-nb_pass+1,-1
                    DO i = 1, nx
                        DO p = 1, npass 

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO 
                    END DO
                END DO
            END DO 
        END IF
        
        CASE(9) ! x-edge 3

        DO k = nz, nz-nb+1, -1  
            DO j = 1, nb        
                DO i = 1, nx
                    DO p = 1, 7

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO 

        DO k = nz, nz-nb+1, -1  
            DO j = 1, nb        
                DO i = 1, nx
                    DO p = 1, 3

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO 

        IF(npass .GT. 0) THEN
            DO k = nz, nz-nb_pass+1, -1  
                DO j = 1, nb_pass        
                    DO i = 1, nx
                        DO p = 1, npass

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO 
                    END DO
                END DO
            END DO 
        END IF
        
        CASE(10) ! x-edge 4

        DO k = nz, nz-nb+1, -1
            DO j = ny, ny-nb+1, -1
                DO i = 1, nx
                    DO p = 1, 7 

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO    
        END DO
        
        DO k = nz, nz-nb+1, -1
            DO j = ny, ny-nb+1, -1
                DO i = 1, nx
                    DO p = 1, 3 

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO    
        END DO

        IF(npass .GT. 0) THEN
            DO k = nz, nz-nb_pass+1, -1
                DO j = ny, ny-nb_pass+1, -1
                    DO i = 1, nx
                        DO p = 1, npass 

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO 
                    END DO
                END DO    
            END DO
        END IF
    
        CASE(11) ! y-edge 1

        DO k = 1, nb
            DO i = 1, nb
                DO j = 1, ny
                    DO p = 1, 7 

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO    

        DO k = 1, nb
            DO i = 1, nb
                DO j = 1, ny
                    DO p = 1, 3 

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO    

        IF(npass .GT. 0) THEN
            DO k = 1, nb_pass
                DO i = 1, nb_pass
                    DO j = 1, ny
                        DO p = 1, npass 

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO 
                    END DO
                END DO
            END DO 
        END IF
        
        CASE(12) ! y-edge 2

        DO k = 1, nb
            DO i = nx, nx-nb+1, -1
                DO j = 1, ny
                    DO p = 1, 7 

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO    

        DO k = 1, nb
            DO i = nx, nx-nb+1, -1
                DO j = 1, ny
                    DO p = 1, 3

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO    

        IF(npass .GT. 0) THEN
            DO k = 1, nb_pass
                DO i = nx, nx-nb_pass+1, -1
                    DO j = 1, ny
                        DO p = 1, npass 

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO 
                    END DO
                END DO
            END DO 
        END IF
        
        CASE(13) ! y-edge 3

        DO k = nz, nz-nb+1, -1
            DO i = 1, nb 
                DO j = 1, ny
                    DO p = 1, 7 

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO    
        
        DO k = nz, nz-nb+1, -1
            DO i = 1, nb 
                DO j = 1, ny
                    DO p = 1, 3 

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO 

        IF(npass .GT. 0) THEN
            DO k = nz, nz-nb_pass+1, -1
                DO i = 1, nb_pass 
                    DO j = 1, ny
                        DO p = 1, npass 

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO 
                    END DO
                END DO
            END DO   
        END IF
        
        CASE(14) ! y-edge 4

        DO k = nz, nz-nb+1, -1
            DO i = nx, nx-nb+1, -1
                DO j = 1, ny
                    DO p = 1, 7 

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
               END DO 
            END DO      
        END DO        
        
        DO k = nz, nz-nb+1, -1
            DO i = nx, nx-nb+1, -1
                DO j = 1, ny
                    DO p = 1, 3 

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
               END DO 
            END DO      
        END DO        

        IF(npass .GT. 0) THEN
            DO k = nz, nz-nb_pass+1, -1
                DO i = nx, nx-nb_pass+1, -1
                    DO j = 1, ny
                        DO p = 1, npass 

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO 
                   END DO 
                END DO      
            END DO  
        END IF
        
        CASE(15) ! z-edge 1

        DO j = 1, nb
            DO i = 1, nb
                DO k = 1, nz
                    DO p = 1, 7 

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO

        DO j = 1, nb
            DO i = 1, nb
                DO k = 1, nz
                    DO p = 1, 3

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO

        IF(npass .GT. 0) THEN
            DO j = 1, nb_pass
                DO i = 1, nb_pass
                    DO k = 1, nz
                        DO p = 1, npass 

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO 
                    END DO
                END DO
            END DO      
        END IF

        CASE(16) ! z-edge 2

        DO j = 1, nb
            DO i = nx, nx-nb+1, -1
                DO k = 1, nz
                    DO p = 1, 7 

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO

        DO j = 1, nb
            DO i = nx, nx-nb+1, -1
                DO k = 1, nz
                    DO p = 1, 3 

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO

        IF(npass .GT. 0) THEN
            DO j = 1, nb_pass
                DO i = nx, nx-nb_pass+1, -1
                    DO k = 1, nz
                        DO p = 1, npass 

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO 
                    END DO
                END DO
            END DO
        END IF
        
        CASE(17) ! z-edge 3

        DO j = ny, ny-nb+1, -1 
            DO i = 1, nb
                DO k = 1, nz
                    DO p = 1, 7 

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO

        DO j = ny, ny-nb+1, -1 
            DO i = 1, nb
                DO k = 1, nz
                    DO p = 1, 3 

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO

        IF(npass .GT. 0) THEN
            DO j = ny, ny-nb_pass+1, -1 
                DO i = 1, nb_pass
                    DO k = 1, nz
                        DO p = 1, npass 

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO 
                    END DO
                END DO
            END DO
        END IF

        CASE(18) ! z-edge 4

        DO j = ny, ny-nb+1, -1
            DO i = nx, nx-nb+1, -1
                DO k = 1, nz
                    DO p = 1, 7 

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO  
            END DO  
        END DO  

        DO j = ny, ny-nb+1, -1
            DO i = nx, nx-nb+1, -1
                DO k = 1, nz
                    DO p = 1, 3 

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO  
            END DO  
        END DO  
        
        IF(npass .GT. 0) THEN
            DO j = ny, ny-nb_pass+1, -1
                DO i = nx, nx-nb_pass+1, -1
                    DO k = 1, nz
                        DO p = 1, npass 

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO 
                    END DO  
                END DO  
            END DO  
        END IF
        
        !********************
        ! Corners
        !********************
        
        CASE(19) ! bottom SW corner

        DO k = 1, nb
            DO j = 1, nb
                DO i = 1, nb 
                    DO p = 1, 7 

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO                

        DO k = 1, nb
            DO j = 1, nb
                DO i = 1, nb 
                    DO p = 1, 3 

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO    

        IF(npass .GT. 0) THEN
            DO k = 1, nb_pass
                DO j = 1, nb_pass
                    DO i = 1, nb_pass 
                        DO p = 1, npass 

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO 
                    END DO
                END DO
            END DO 
        END IF
        
        CASE(20) ! bottom SE corner

        DO k = 1, nb
            DO j = 1, nb
                DO i = nx, nx-nb+1, -1
                    DO p = 1, 7 

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO

        DO k = 1, nb
            DO j = 1, nb
                DO i = nx, nx-nb+1, -1
                    DO p = 1, 3 

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO

        IF(npass .GT. 0) THEN
            DO k = 1, nb_pass
                DO j = 1, nb_pass
                    DO i = nx, nx-nb_pass+1, -1
                        DO p = 1, npass 

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO 
                    END DO
                END DO
            END DO
        END IF
        
        CASE(21) ! bottom NW corner

        DO k = 1, nb
            DO j = ny, ny-nb+1, -1
                DO i = 1, nb
                    DO p = 1, 7 

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO
        
        DO k = 1, nb
            DO j = ny, ny-nb+1, -1
                DO i = 1, nb
                    DO p = 1, 3 

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO

        IF(npass .GT. 0) THEN
            DO k = 1, nb_pass
                DO j = ny, ny-nb_pass+1, -1
                    DO i = 1, nb_pass
                        DO p = 1, npass 

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO 
                    END DO
                END DO
            END DO
        END IF
        
        CASE(22) ! bottom NE corner

        DO k = 1, nb
            DO j = ny, ny-nb+1, -1
                DO i = nx, nx-nb+1, -1
                    DO p = 1, 7 

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO
        
        DO k = 1, nb
            DO j = ny, ny-nb+1, -1
                DO i = nx, nx-nb+1, -1
                    DO p = 1, 3 

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO

        IF(npass .GT. 0) THEN
            DO k = 1, nb_pass
                DO j = ny, ny-nb_pass+1, -1
                    DO i = nx, nx-nb_pass+1, -1
                        DO p = 1, npass 

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO 
                    END DO
                END DO
            END DO
        END IF
                    
        CASE(23) ! top SW corner

        DO k = nz, nz-nb+1, -1
            DO j = 1, nb
                DO i = 1, nb
                    DO p = 1, 7 

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO
           
        DO k = nz, nz-nb+1, -1
            DO j = 1, nb
                DO i = 1, nb
                    DO p = 1, 3

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO

        IF(npass .GT. 0) THEN
            DO k = nz, nz-nb_pass+1, -1
                DO j = 1, nb_pass
                    DO i = 1, nb_pass
                        DO p = 1, npass 

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO 
                    END DO
                END DO
            END DO
        END IF
        
        CASE(24) ! top SE corner
    
        DO k = nz, nz-nb+1, -1
            DO j = 1, nb
                DO i = nx, nx-nb+1, -1
                    DO p = 1, 7 

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO

        DO k = nz, nz-nb+1, -1
            DO j = 1, nb
                DO i = nx, nx-nb+1, -1
                    DO p = 1, 3 

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO

        IF(npass .GT. 0) THEN
            DO k = nz, nz-nb_pass+1, -1
                DO j = 1, nb_pass
                    DO i = nx, nx-nb_pass+1, -1
                        DO p = 1, npass

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO 
                    END DO
                END DO
            END DO
        END IF

        CASE(25) ! top NW corner

        DO k = nz, nz-nb+1, -1
            DO j = ny, ny-nb+1, -1
                DO i = 1, nb
                    DO p = 1, 7 

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO
            
        DO k = nz, nz-nb+1, -1
            DO j = ny, ny-nb+1, -1
                DO i = 1, nb
                    DO p = 1, 3 

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO

        IF(npass .GT. 0) THEN
            DO k = nz, nz-nb_pass+1, -1
                DO j = ny, ny-nb_pass+1, -1
                    DO i = 1, nb_pass
                        DO p = 1, npass 

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO 
                    END DO
                END DO
            END DO
        END IF

        CASE(26) ! top NE corner

        DO k = nz, nz-nb+1, -1
            DO j = ny, ny-nb+1, -1
                DO i = nx, nx-nb+1, -1
                    DO p = 1, 7 

                        buff(ix) = q_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO
         
        DO k = nz, nz-nb+1, -1
            DO j = ny, ny-nb+1, -1
                DO i = nx, nx-nb+1, -1
                    DO p = 1, 3 

                        buff(ix) = bface_3d(i,j,k,p)
                        ix = ix +1
 
                    END DO 
                END DO
            END DO
        END DO

        IF(npass .GT. 0) THEN
            DO k = nz, nz-nb_pass+1, -1
                DO j = ny, ny-nb_pass+1, -1
                    DO i = nx, nx-nb_pass+1, -1
                        DO p = 1, npass

                            buff(ix) = qpass_3d(i,j,k,p)
                            ix = ix +1
     
                        END DO 
                    END DO
                END DO
            END DO
        END IF
        
     END SELECT   

END SUBROUTINE packup_outbuffer


SUBROUTINE unpack_bufferin()

    INTEGER :: i, j, k, p, ix, nb_pass

    nb_pass = 3

    !********************
    ! Faces
    !********************
         
    ! x- boundary
    
    IF((neighbor_rank(1) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(1) .NE. myrank)) THEN

    ix = 2 + 1 + disp(1)
    
    DO i = 0, 1-nb, -1
        DO k = 1, nz
            DO j = 1, ny
                DO p = 1, 7         
                                              
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1    
                
                END DO
            END DO
        END DO
    END DO
    
    DO i = 0, 1-nb, -1
        DO k = 1, nz
            DO j = 1, ny
                DO p = 1, 3         
                                              
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1    
                
                END DO
            END DO
        END DO
    END DO
            
    IF(npass .GT. 0) THEN
        DO i = 0, 1-nb_pass, -1
            DO k = 1, nz
                DO j = 1, ny
                    DO p = 1, npass         
                                                  
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1    
                    
                    END DO
                END DO
            END DO
        END DO
    END IF
            
    END IF

    IF((neighbor_rank(2) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(2) .NE. myrank)) THEN
            
    ! x+ boundary
    ix = 2 + 1 + disp(2)
    
    DO i = nx+1, nx+nb
        DO k = 1, nz
            DO j = 1, ny
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO
   
    DO i = nx+1, nx+nb
        DO k = 1, nz
            DO j = 1, ny
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO

    IF(npass .GT. 0) THEN
        DO i = nx+1, nx+nb_pass
            DO k = 1, nz
                DO j = 1, ny
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO
            END DO
        END DO
    END IF
    
    END IF
    
    IF((neighbor_rank(3) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(3) .NE. myrank)) THEN

    ! y- boundary
    ix = 2 + 1 + disp(3)
    
    DO j = 0, 1-nb, -1 
        DO k = 1, nz
            DO i = 1, nx
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO        
        END DO
    END DO

    DO j = 0, 1-nb, -1 
        DO k = 1, nz
            DO i = 1, nx
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO        
        END DO
    END DO

    IF(npass .GT. 0) THEN
        DO j = 0, 1-nb_pass, -1 
            DO k = 1, nz
                DO i = 1, nx
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO        
            END DO
        END DO
    END IF
    
    END IF

    IF((neighbor_rank(4) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(4) .NE. myrank)) THEN
    
    ! y+ boundary
    ix = 2 + 1 + disp(4)
    
    DO j = ny+1, ny+nb
        DO k = 1, nz
            DO i = 1, nx
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO            
        END DO
    END DO       

    DO j = ny+1, ny+nb
        DO k = 1, nz
            DO i = 1, nx
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO            
        END DO
    END DO       

    IF(npass .GT. 0) THEN
        DO j = ny+1, ny+nb_pass
            DO k = 1, nz
                DO i = 1, nx
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO            
            END DO
        END DO  
    END IF
    
    END IF
    
    IF((neighbor_rank(5) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(5) .NE. myrank)) THEN

    ! z- boundary
    ix = 2 + 1 + disp(5)
    
    DO k = 0, 1-nb, -1
        DO j = 1, ny
            DO i = 1, nx
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO    
        END DO
    END DO
    
    DO k = 0, 1-nb, -1
        DO j = 1, ny
            DO i = 1, nx
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO    
        END DO
    END DO

    IF(npass .GT. 0) THEN
        DO k = 0, 1-nb_pass, -1
            DO j = 1, ny
                DO i = 1, nx
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO    
            END DO
        END DO
    END IF
    
    END IF

    IF((neighbor_rank(6) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(6) .NE. myrank)) THEN
    
    ! z+ boundary
    ix = 2 + 1 + disp(6)
    
    DO k = nz+1, nz+nb
        DO j = 1, ny
            DO i = 1, nx
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO    
        END DO
    END DO
    
    DO k = nz+1, nz+nb
        DO j = 1, ny
            DO i = 1, nx
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO    
        END DO
    END DO

    IF(npass .GT. 0) THEN
        DO k = nz+1, nz+nb_pass
            DO j = 1, ny
                DO i = 1, nx
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO    
            END DO
        END DO
    END IF
    
    END IF
    
    !********************
    ! Edges
    !********************
    
    IF((neighbor_rank(7) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(7) .NE. myrank)) THEN

    ! x-edge 1
    ix = 2 + 1 + disp(7)
    
    DO k = 0, 1-nb, -1
        DO j = 0, 1-nb, -1
            DO i = 1, nx
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO    
        END DO            
    END DO

    DO k = 0, 1-nb, -1
        DO j = 0, 1-nb, -1
            DO i = 1, nx
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO    
        END DO            
    END DO

    IF(npass .GT. 0) THEN
        DO k = 0, 1-nb_pass, -1
            DO j = 0, 1-nb_pass, -1
                DO i = 1, nx
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO    
            END DO            
        END DO
    END IF
    
    END IF
    
    IF((neighbor_rank(8) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(8) .NE. myrank)) THEN

    ! x-edge 2
    ix = 2 + 1 + disp(8)
    
    DO k = 0, 1-nb, -1
        DO j = ny+1, ny+nb
            DO i = 1, nx
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO    
        END DO        
    END DO
 
    DO k = 0, 1-nb, -1
        DO j = ny+1, ny+nb
            DO i = 1, nx
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO    
        END DO        
    END DO
    
    IF(npass .GT. 0) THEN
        DO k = 0, 1-nb_pass, -1
            DO j = ny+1, ny+nb_pass
                DO i = 1, nx
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO    
            END DO        
        END DO
    END IF    
    
    END IF
   
    IF((neighbor_rank(9) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(9) .NE. myrank)) THEN

    ! x-edge 3
    ix = 2 + 1 + disp(9)
    
    DO k = nz+1, nz+nb
        DO j = 0, 1-nb, -1
            DO i = 1, nx
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO    
        END DO        
    END DO

    DO k = nz+1, nz+nb
        DO j = 0, 1-nb, -1
            DO i = 1, nx
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO    
        END DO        
    END DO

    IF(npass .GT. 0) THEN
        DO k = nz+1, nz+nb_pass
            DO j = 0, 1-nb_pass, -1
                DO i = 1, nx
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO    
            END DO        
        END DO
    END IF
    
    END IF
    
    IF((neighbor_rank(10) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(10) .NE. myrank)) THEN

    ! x-edge 4
    ix = 2 + 1 + disp(10)
    
    DO k = nz+1, nz+nb
        DO j = ny+1, ny+nb
            DO i = 1, nx
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO            
    END DO

    DO k = nz+1, nz+nb
        DO j = ny+1, ny+nb
            DO i = 1, nx
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO            
    END DO

    IF(npass .GT. 0) THEN
        DO k = nz+1, nz+nb_pass
            DO j = ny+1, ny+nb_pass
                DO i = 1, nx
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO
            END DO            
        END DO
    END IF
    
    END IF
    
    IF((neighbor_rank(11) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(11) .NE. myrank)) THEN

    ! y-edge 1
    ix = 2 + 1 + disp(11)
    
    DO k = 0, 1-nb, -1
        DO i = 0, 1-nb, -1
            DO j = 1, ny
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO    
        END DO    
    END DO

    DO k = 0, 1-nb, -1
        DO i = 0, 1-nb, -1
            DO j = 1, ny
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO    
        END DO    
    END DO

    IF(npass .GT. 0) THEN
        DO k = 0, 1-nb_pass, -1
            DO i = 0, 1-nb_pass, -1
                DO j = 1, ny
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO    
            END DO    
        END DO
    END IF
    
    END IF

    IF((neighbor_rank(12) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(12) .NE. myrank)) THEN

    ! y-edge 2
    ix = 2 + 1 + disp(12)
    
    DO k = 0, 1-nb, -1
        DO i = nx+1, nx+nb
            DO j = 1, ny
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO
    
    DO k = 0, 1-nb, -1
        DO i = nx+1, nx+nb
            DO j = 1, ny
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO
    
    IF(npass .GT. 0) THEN
        DO k = 0, 1-nb_pass, -1
            DO i = nx+1, nx+nb_pass
                DO j = 1, ny
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO
            END DO
        END DO
    END IF
    
    END IF
    
    IF((neighbor_rank(13) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(13) .NE. myrank)) THEN

    ! y-edge 3
    ix = 2 + 1 + disp(13)
    
    DO k =  nz+1, nz+nb
        DO i = 0, 1-nb ,-1 
            DO j = 1, ny
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO
   
    DO k =  nz+1, nz+nb
        DO i = 0, 1-nb ,-1 
            DO j = 1, ny
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO

    IF(npass .GT. 0) THEN
        DO k = nz+1, nz+nb_pass
            DO i = 0, 1-nb_pass ,-1 
                DO j = 1, ny
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO
            END DO
        END DO
    END IF
    
    END IF
    
    IF((neighbor_rank(14) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(14) .NE. myrank)) THEN

    ! y-edge 4
    ix = 2 + 1 + disp(14)
    
    DO k = nz+1, nz+nb
        DO i = nx+1, nx+nb
            DO j = 1, ny
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO
   
    DO k = nz+1, nz+nb
        DO i = nx+1, nx+nb
            DO j = 1, ny
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO

    IF(npass .GT. 0) THEN
        DO k = nz+1, nz+nb_pass
            DO i = nx+1, nx+nb_pass
                DO j = 1, ny
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO
            END DO
        END DO
    END IF
    
    END IF

    IF((neighbor_rank(15) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(15) .NE. myrank)) THEN

    ! z-edge 1
    ix = 2 + 1 + disp(15)
    
    DO j = 0, 1-nb, -1
        DO i = 0, 1-nb, -1 
            DO k = 1, nz
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO
    
    DO j = 0, 1-nb, -1
        DO i = 0, 1-nb, -1 
            DO k = 1, nz
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO

    IF(npass .GT. 0) THEN
        DO j = 0, 1-nb_pass, -1
            DO i = 0, 1-nb_pass, -1 
                DO k = 1, nz
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO
            END DO
        END DO
    END IF
    
    END IF

    IF((neighbor_rank(16) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(16) .NE. myrank)) THEN

    ! z-edge 2
    ix = 2 + 1 + disp(16)
    
    DO j = 0, 1-nb, -1
        DO i = nx+1, nx+nb
            DO k = 1, nz
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO
    
    DO j = 0, 1-nb, -1
        DO i = nx+1, nx+nb
            DO k = 1, nz
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO

    IF(npass .GT. 0) THEN
        DO j = 0, 1-nb_pass, -1
            DO i = nx+1, nx+nb_pass
                DO k = 1, nz
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO
            END DO
        END DO
    END IF
    
    END IF

    IF((neighbor_rank(17) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(17) .NE. myrank)) THEN

    ! z-edge 3
    ix = 2 + 1 + disp(17)
    
    DO j = ny+1, ny+nb
        DO i = 0, 1-nb, -1
            DO k = 1, nz
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO
    
    DO j = ny+1, ny+nb
        DO i = 0, 1-nb, -1
            DO k = 1, nz
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO

    IF(npass .GT. 0) THEN
        DO j = ny+1, ny+nb_pass
            DO i = 0, 1-nb_pass, -1
                DO k = 1, nz
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO
            END DO
        END DO
    END IF
    
    END IF
    
    IF((neighbor_rank(18) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(18) .NE. myrank)) THEN

    ! z-edge 4
    ix = 2 + 1 + disp(18)
    
    DO j = ny+1, ny+nb
        DO i = nx+1, nx+nb
            DO k = 1, nz
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO
    
    DO j = ny+1, ny+nb
        DO i = nx+1, nx+nb
            DO k = 1, nz
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO

    IF(npass .GT. 0) THEN
        DO j = ny+1, ny+nb_pass
            DO i = nx+1, nx+nb_pass
                DO k = 1, nz
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO
            END DO
        END DO
    END IF
    
    END IF
    
    !********************
    ! Corners
    !********************
    
    IF((neighbor_rank(19) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(19) .NE. myrank)) THEN

    ! bottom SW corner 
    ix = 2 + 1 + disp(19)
    
    DO k = 0, 1-nb, -1
        DO j = 0, 1-nb, -1
            DO i = 0, 1-nb, -1
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO

    DO k = 0, 1-nb, -1
        DO j = 0, 1-nb, -1
            DO i = 0, 1-nb, -1
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO

    IF(npass .GT. 0) THEN
        DO k = 0, 1-nb_pass, -1
            DO j = 0, 1-nb_pass, -1
                DO i = 0, 1-nb_pass, -1
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO
            END DO
        END DO
    END IF
    
    END IF
    
    IF((neighbor_rank(20) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(20) .NE. myrank)) THEN

    ! bottom SE corner 
    ix = 2 + 1 + disp(20)

    DO k = 0, 1-nb, -1
        DO j = 0, 1-nb, -1
            DO i = nx+1, nx+nb
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO

    DO k = 0, 1-nb, -1
        DO j = 0, 1-nb, -1
            DO i = nx+1, nx+nb
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO

    IF(npass .GT. 0) THEN
        DO k = 0, 1-nb_pass, -1
            DO j = 0, 1-nb_pass, -1
                DO i = nx+1, nx+nb_pass
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO
            END DO
        END DO
    END IF
    
    END IF

    IF((neighbor_rank(21) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(21) .NE. myrank)) THEN

    ! bottom NW corner 
    ix = 2 + 1 + disp(21)
    
    DO k = 0, 1-nb, -1
        DO j = ny+1, ny+nb
            DO i = 0, 1-nb, -1
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO
    
    DO k = 0, 1-nb, -1
        DO j = ny+1, ny+nb
            DO i = 0, 1-nb, -1
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO

    IF(npass .GT. 0) THEN

    END IF
        DO k = 0, 1-nb_pass, -1
            DO j = ny+1, ny+nb_pass
                DO i = 0, 1-nb_pass, -1
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO
            END DO
        END DO
    END IF

    IF((neighbor_rank(22) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(22) .NE. myrank)) THEN

    ! bottom NE corner 
    ix = 2 + 1 + disp(22)

    DO k = 0, 1-nb, -1
        DO j = ny+1, ny+nb
            DO i = nx+1, nx+nb
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO

    DO k = 0, 1-nb, -1
        DO j = ny+1, ny+nb
            DO i = nx+1, nx+nb
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO

    IF(npass .GT. 0) THEN
        DO k = 0, 1-nb_pass, -1
            DO j = ny+1, ny+nb_pass
                DO i = nx+1, nx+nb_pass
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO
            END DO
        END DO
    END IF
    
    END IF
    
    IF((neighbor_rank(23) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(23) .NE. myrank)) THEN

    ! top SW corner 
    ix = 2 + 1 + disp(23)
    
    DO k = nz+1, nz+nb
        DO j = 0, 1-nb, -1
            DO i = 0, 1-nb,-1
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO

    DO k = nz+1, nz+nb
        DO j = 0, 1-nb, -1
            DO i = 0, 1-nb,-1
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO

    IF(npass .GT. 0) THEN
        DO k = nz+1, nz+nb_pass
            DO j = 0, 1-nb_pass, -1
                DO i = 0, 1-nb_pass,-1
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO
            END DO
        END DO
    END IF
    
    END IF
    
    IF((neighbor_rank(24) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(24) .NE. myrank)) THEN

    ! top SE corner 
    ix = 2 + 1 + disp(24)
    
    DO k = nz+1, nz+nb
        DO j = 0, 1-nb, -1
            DO i = nx+1, nx+nb
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO

    DO k = nz+1, nz+nb
        DO j = 0, 1-nb, -1
            DO i = nx+1, nx+nb
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO

    IF(npass .GT. 0) THEN
        DO k = nz+1, nz+nb_pass
            DO j = 0, 1-nb_pass, -1
                DO i = nx+1, nx+nb_pass
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO
            END DO
        END DO
    END IF
    
    END IF
    
    IF((neighbor_rank(25) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(25) .NE. myrank)) THEN

    ! top NW corner 
    ix = 2 + 1 + disp(25)
    
    DO k = nz+1, nz+nb
        DO j = ny+1, ny+nb
            DO i = 0, 1-nb, -1
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO

    DO k = nz+1, nz+nb
        DO j = ny+1, ny+nb
            DO i = 0, 1-nb, -1
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO

    IF(npass .GT. 0) THEN
        DO k = nz+1, nz+nb_pass
            DO j = ny+1, ny+nb_pass
                DO i = 0, 1-nb_pass, -1
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO
            END DO
        END DO
    END IF
    
    END IF
    
    IF((neighbor_rank(26) .NE. MPI_PROC_NULL) .AND. (neighbor_rank(26) .NE. myrank)) THEN

    ! top NE corner 
    ix = 2 + 1 + disp(26)
    
    DO k = nz+1, nz+nb
        DO j = ny+1, ny+nb
            DO i = nx+1, nx+nb
                DO p = 1, 7
          
                    q_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO
    
    DO k = nz+1, nz+nb
        DO j = ny+1, ny+nb
            DO i = nx+1, nx+nb
                DO p = 1, 3
          
                    bface_3d(i,j,k,p) = A(ix) 
                    ix = ix +1
                
                END DO
            END DO
        END DO
    END DO

    IF(npass .GT. 0) THEN
        DO k = nz+1, nz+nb_pass
            DO j = ny+1, ny+nb_pass
                DO i = nx+1, nx+nb_pass
                    DO p = 1, npass
              
                        qpass_3d(i,j,k,p) = A(ix) 
                        ix = ix +1
                    
                    END DO
                END DO
            END DO
        END DO
    END IF
    
    END IF
        
    
END SUBROUTINE unpack_bufferin


END MODULE mpi_domain_mod