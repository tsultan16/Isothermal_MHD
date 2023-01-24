! This example program demonstrates halo exchange with 3d domain decomposition
! using one-sided RMA. Each process directly puts data from it's top and 
! bottom interior layers into the neighbor's RMA window. 
! We use passive target RMA and acheive synchronization via MPI_WIN_FLUSH. 


PROGRAM test7_onesided

USE ISO_C_BINDING  ! needed for access to C pointer
!USE MPI

IMPLICIT NONE

INCLUDE 'mpif.h'  ! MPI runtime library access


!***************************************************************************************

INTEGER, PARAMETER :: nx = 3
INTEGER, PARAMETER :: ny = nx, nz = nx  
INTEGER, PARAMETER :: nb = 1
INTEGER, PARAMETER :: ndim = 3
INTEGER, PARAMETER :: nranks_x = 2
INTEGER, PARAMETER :: nranks_y = 2
INTEGER, PARAMETER :: nranks_z = 2

INTEGER :: comm3d, myrank, dims(3), numprocs(1), coord(3), req(6), status(MPI_STATUS_SIZE)
LOGICAL :: isperiodic(3), reorder 
INTEGER :: xlow, xhi, ylow, yhi, zlow, zhi, neighbor, ierr, neighbor_rank(26)
REAL(8), ALLOCATABLE :: u(:,:,:), buffer(:)
INTEGER :: i, j, k, np

!RMA Vars
INTEGER :: win, sizedouble, winsize
INTEGER(KIND=MPI_ADDRESS_KIND) :: disp(26) ! 6 faces + 12 edges + 8 corners
INTEGER(KIND=MPI_ADDRESS_KIND) :: shuffled_disp(26)
INTEGER :: send_buffer_size(26)
INTEGER :: window
TYPE(C_PTR) :: rma_cmem  ! C pointer to RMA window memory block
REAL(8), POINTER :: A(:) ! Fortran pointer that will be associated with the RMA window C pointer

!***************************************************************************************


! allocate work arrays
ALLOCATE(u(1-nb:nx+nb, 1-nb:ny+nb, 1-nb:nz+nb), buffer(MAX(nx*ny,ny*nz,nz*nx)))
u = 0.d0
buffer = 0.d0


! MPI initialization
CALL MPI_INIT(ierr)

CALL mpi_comm3d_init()

CALL mpi_get_neighbor_ranks()

CALL rma_init()

CALL MPI_BARRIER(comm3d, ierr)

!GO TO 111

! open passive RMA exposure epoch
CALL open_exposure_epoch()

CALL MPI_BARRIER(comm3d, ierr)

! put some test values inside u
u =  myrank

CALL print_array(0)

! tell neighbors that we are ready to recv
CALL recv_signal()

! poll on signal neighbors telling us they are ready to recv
CALL poll_signal(A,1)

! send data to neighbors
CALL send_msg()

! tell neighbors that we are done sending
CALL data_sent_signal()

! poll on signal from from neighbors telling us they have sent us data_sent_signal
CALL poll_signal(A,2)


PRINT*,''
PRINT*,'One sided communication completed...myrank = ',myrank
PRINT*,''


CALL MPI_BARRIER(comm3d, ierr)


!unpack data from recv buffer
CALL unpack_bufferin()


CALL print_array(1)



CALL close_rma_window()


!111 CONTINUE


! shut off mpi
CALL MPI_FINALIZE(ierr)

PRINT*,'Done..myrank = ', myrank


DEALLOCATE(u, buffer)


STOP


CONTAINS



! send signal to neighbor indicating that we are ready to receive data
SUBROUTINE recv_signal()

    INTEGER :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND) :: msg_disp

    buffer(1) = 1.d0
    
    ! Put recv signal in neighbor rma windows                 
    DO i = 1, 26             
   
        msg_disp = disp(i)    
        CALL MPI_PUT(buffer, 1, MPI_DOUBLE_PRECISION, neighbor_rank(i), &
                 msg_disp, 1, MPI_DOUBLE_PRECISION, win, ierr)            
         
    END DO  
    
    CALL MPI_WIN_FLUSH(neighbor, win, ierr)


END SUBROUTINE recv_signal


! send a signal to neighbor indicating that we have sent them data
SUBROUTINE data_sent_signal()

    INTEGER :: ierr, i
    INTEGER(KIND=MPI_ADDRESS_KIND) :: msg_disp

    buffer(1) = 1.d0

    ! Put sent signal in neighbor rma windows                 
    DO i = 1, 26             
   
        msg_disp = 1 + disp(i)    
        CALL MPI_PUT(buffer, 1, MPI_DOUBLE_PRECISION, neighbor_rank(i), &
                 msg_disp, 1, MPI_DOUBLE_PRECISION, win, ierr)            
         
    END DO             
                 
    CALL MPI_WIN_FLUSH(neighbor, win, ierr)


END SUBROUTINE data_sent_signal


! polls on signal from neighbors
SUBROUTINE poll_signal(rma_window, signal_offset)

    REAL(8), VOLATILE :: rma_window(:) ! must include the volatile attribute here.. to ensure that compiler isn't using a copy of the RMA window that it has pre-loaded into a register
    INTEGER, INTENT(IN) :: signal_offset  ! 1 for recv signal, 2 for sent signal
    INTEGER :: ierr, i, poll_stat
    LOGICAL :: poll_flag
    
    poll_flag = .FALSE.
    
    ! poll the signal status
    DO WHILE(.NOT. poll_flag)
    
        poll_flag = .TRUE.
   
        DO i = 1,26
             poll_flag = poll_flag .AND. (rma_window(signal_offset + disp(i)) .EQ. 1.d0)
        END DO
    
        ! wait till signals from all neighbors have arrived...
    
    END DO

    ! reset the signal status
    DO i = 1,26
        rma_window(signal_offset + disp(i)) = 0.d0 
    END DO

END SUBROUTINE poll_signal


SUBROUTINE send_msg()

    INTEGER :: ierr, i, buff_size
    INTEGER(KIND=MPI_ADDRESS_KIND) :: msg_disp


    ! put boundary data in neighbor window
    DO i = 1, 26
    
        msg_disp = 2 + shuffled_disp(i)    
        buff_size = send_buffer_size(i)        
        
        CALL packup_outbuffer(i) 
        
        CALL MPI_PUT(buffer, buff_size, MPI_DOUBLE_PRECISION, neighbor_rank(i), &
                 msg_disp, buff_size, MPI_DOUBLE_PRECISION, win, ierr)
               
    END DO

    ! force all PUTs to complete           
    CALL MPI_WIN_FLUSH(neighbor, win, ierr)


END SUBROUTINE send_msg


SUBROUTINE packup_outbuffer(bndry)

    INTEGER, INTENT(IN) :: bndry
    INTEGER :: i, j, k, ix

    ix = 1

     SELECT CASE(bndry)

        CASE(1) ! x- boundary
          
            DO k = 1, nz
            DO j = 1, ny
          
                buffer(ix) = u(1,j,k)
                ix = ix +1
                
            END DO
            END DO
            
        CASE(2) ! x+ boundary

            DO k = 1, nz
            DO j = 1, ny
          
                buffer(ix) = u(nx,j,k)
                ix = ix +1
                
            END DO
            END DO
        
        CASE(3) ! y- boundary

            DO k = 1, nz
            DO i = 1, nx
          
                buffer(ix) = u(i,1,k)
                ix = ix +1
                
            END DO
            END DO
        
        CASE(4) ! y+ boundary

            DO k = 1, nz
            DO i = 1, nx
          
                buffer(ix) = u(i,ny,k)
                ix = ix +1
                
            END DO
            END DO
      
        CASE(5) ! z- boundary

            DO j = 1, ny
            DO i = 1, nx
          
                buffer(ix) = u(i,j,1)
                ix = ix +1
                
            END DO
            END DO
        
        CASE(6) ! z+ boundary

            DO j = 1, ny
            DO i = 1, nx
          
                buffer(ix) = u(i,j,nz)
                ix = ix +1
                
            END DO
            END DO
        

        CASE(7) ! x-edge 1

            DO i = 1, nx
          
            buffer(ix) = u(i,1,1)
            ix = ix +1
                
            END DO
        
        CASE(8) ! x-edge 2

            DO i = 1, nx
          
            buffer(ix) = u(i,ny,1)
            ix = ix +1
                
            END DO
        
        CASE(9) ! x-edge 3

            DO i = 1, nx
          
            buffer(ix) = u(i,1,nz)
            ix = ix +1
                
            END DO
        
        CASE(10) ! x-edge 4

            DO i = 1, nx
          
            buffer(ix) = u(i,ny,nz)
            ix = ix +1
                
            END DO
        
        CASE(11) ! y-edge 1

            DO j = 1, ny
          
            buffer(ix) = u(1,j,1)
            ix = ix +1
                
            END DO
        
        CASE(12) ! y-edge 2

            DO j = 1, ny
          
            buffer(ix) = u(nx,j,1)
            ix = ix +1
                
            END DO
     
        CASE(13) ! y-edge 3

            DO j = 1, ny
          
            buffer(ix) = u(1,j,nz)
            ix = ix +1
                
            END DO
        
        CASE(14) ! y-edge 4

            DO j = 1, ny
          
            buffer(ix) = u(nx,j,nz)
            ix = ix +1
                
            END DO     
        
         CASE(15) ! z-edge 1

            DO k = 1, nz
          
            buffer(ix) = u(1,1,k)
            ix = ix +1
                
            END DO

         CASE(16) ! z-edge 2

            DO k = 1, nz
          
            buffer(ix) = u(nx,1,k)
            ix = ix +1
                
            END DO
          
         CASE(17) ! z-edge 3

            DO k = 1, nz
          
            buffer(ix) = u(1,ny,k)
            ix = ix +1
                
            END DO

         CASE(18) ! z-edge 4

            DO k = 1, nz
          
            buffer(ix) = u(nx,ny,k)
            ix = ix +1
                
            END DO  

         CASE(19) ! bottom SW corner

            buffer(ix) = u(1,1,1)
            ix = ix +1
                
         CASE(20) ! bottom SE corner

            buffer(ix) = u(nx,1,1)
            ix = ix +1

         CASE(21) ! bottom NW corner

            buffer(ix) = u(1,ny,1)
            ix = ix +1

         CASE(22) ! bottom NE corner

            buffer(ix) = u(nx,ny,1)
            ix = ix +1
            
        CASE(23) ! top SW corner

            buffer(ix) = u(1,1,nz)
            ix = ix +1
                
         CASE(24) ! top SE corner

            buffer(ix) = u(nx,1,nz)
            ix = ix +1

         CASE(25) ! top NW corner

            buffer(ix) = u(1,ny,nz)
            ix = ix +1

         CASE(26) ! top NE corner

            buffer(ix) = u(nx,ny,nz)
            ix = ix +1
         
         
     END SELECT   

END SUBROUTINE packup_outbuffer


SUBROUTINE unpack_bufferin()


    INTEGER :: i, j, k, ix

    ! x- boundary
    ix = 2 + 1 + disp(1)
    
    DO k = 1, nz
        DO j = 1, ny
          
           u(1-nb,j,k) = A(ix) 
           ix = ix +1
                
        END DO
    END DO
            
    ! x+ boundary
    ix = 2 + 1 + disp(2)
    
    DO k = 1, nz
        DO j = 1, ny
          
            u(nx+nb,j,k) = A(ix) 
            ix = ix +1    
        
        END DO
    END DO
    
    ! y- boundary
    ix = 2 + 1 + disp(3)
    
    DO k = 1, nz
        DO i = 1, nx
          
            u(i,1-nb,k) = A(ix) 
            ix = ix +1
                
        END DO
    END DO
        
    ! y+ boundary
    ix = 2 + 1 + disp(4)
    
    DO k = 1, nz
        DO i = 1, nx
          
            u(i,ny+nb,k) = A(ix) 
            ix = ix +1
                
        END DO
    END DO       

    ! z- boundary
    ix = 2 + 1 + disp(5)
    
    DO j = 1, ny
        DO i = 1, nx
          
            u(i,j,1-nb) = A(ix) 
            ix = ix +1
                
        END DO
    END DO
        
    ! z+ boundary
    ix = 2 + 1 + disp(6)
    
    DO j = 1, ny
        DO i = 1, nx
          
            u(i,j,nz+nb) = A(ix) 
            ix = ix +1
                
        END DO
    END DO
    
    ! x-edge 1
    ix = 2 + 1 + disp(7)
    
    DO i = 1, nx
          
      u(i,1-nb,1-nb) = A(ix)
      ix = ix +1
                
    END DO
    
    ! x-edge 2
    ix = 2 + 1 + disp(8)
    
    DO i = 1, nx
          
      u(i,ny+nb,1-nb) = A(ix)
      ix = ix +1
                
    END DO
   
    ! x-edge 3
    ix = 2 + 1 + disp(9)
    
    DO i = 1, nx
          
      u(i,1-nb,nz+nb) = A(ix)
      ix = ix +1
                
    END DO
    
    ! x-edge 4
    ix = 2 + 1 + disp(10)
    
    DO i = 1, nx
          
      u(i,ny+nb,nz+nb) = A(ix)
      ix = ix +1
                
    END DO
    
    ! y-edge 1
    ix = 2 + 1 + disp(11)
    
    DO j = 1, ny
          
      u(1-nb,j,1-nb) = A(ix)
      ix = ix +1
                
    END DO
    
    ! y-edge 2
    ix = 2 + 1 + disp(12)
    
    DO j = 1, ny
          
      u(nx+nb,j,1-nb) = A(ix)
      ix = ix +1
                
    END DO
    
    ! y-edge 3
    ix = 2 + 1 + disp(13)
    
    DO j = 1, ny
          
      u(1-nb,j,nz+nb) = A(ix)
      ix = ix +1
                
    END DO
    
    ! y-edge 4
    ix = 2 + 1 + disp(14)
    
    DO j = 1, ny
          
      u(nx+nb,j,nz+nb) = A(ix)
      ix = ix +1
                
    END DO
    
    ! z-edge 1
    ix = 2 + 1 + disp(15)
    
    DO k = 1, nz
          
      u(1-nb,1-nb,k) = A(ix)
      ix = ix +1
                
    END DO
    
    ! z-edge 2
    ix = 2 + 1 + disp(16)
    
    DO k = 1, nz
          
      u(nx+nb,1-nb,k) = A(ix)
      ix = ix +1
                
    END DO
    
    ! z-edge 3
    ix = 2 + 1 + disp(17)
    
    DO k = 1, nz
          
      u(1-nb,ny+nb,k) = A(ix)
      ix = ix +1
                
    END DO
    
    ! z-edge 4
    ix = 2 + 1 + disp(18)
    
    DO k = 1, nz
          
      u(nx+nb,ny+nb,k) = A(ix)
      ix = ix +1
                
    END DO
    
    ! bottom SW corner 
    ix = 2 + 1 + disp(19)
    
    u(1-nb,1-nb,1-nb) = A(ix)

    ! bottom SE corner 
    ix = 2 + 1 + disp(20)
    
    u(nx+nb,1-nb,1-nb) = A(ix)

    ! bottom NW corner 
    ix = 2 + 1 + disp(21)
    
    u(1-nb,ny+nb,1-nb) = A(ix)

    ! bottom NE corner 
    ix = 2 + 1 + disp(22)
    
    u(nx+nb,ny+nb,1-nb) = A(ix)

    ! top SW corner 
    ix = 2 + 1 + disp(23)
    
    u(1-nb,1-nb,nz+nb) = A(ix)

    ! top SE corner 
    ix = 2 + 1 + disp(24)
    
    u(nx+nb,1-nb,nz+nb) = A(ix)

    ! top NW corner 
    ix = 2 + 1 + disp(25)
    
    u(1-nb,ny+nb,nz+nb) = A(ix)

    ! top NE corner 
    ix = 2 + 1 + disp(26)
    
    u(nx+nb,ny+nb,nz+nb) = A(ix)
    
    
    
END SUBROUTINE unpack_bufferin


SUBROUTINE mpi_comm3d_init()

    INTEGER :: ierr
    
    ! find out number of processors
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs(1), ierr)
    
    np = numprocs(1)
    
    !PRINT*,'# of processors =', np
    
    ! create a cartesian communicator for our 3d domain decomposition
    dims(1) = nranks_x
    dims(2) = nranks_y
    dims(3) = nranks_z
    isperiodic(1:3) = .TRUE.  ! periodic boundaries
    reorder = .TRUE.
    
    CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndim, dims, isperiodic, reorder, comm3d, ierr)
    
    ! get rank of this processor
    CALL MPI_COMM_RANK(comm3d, myrank, ierr)

    IF(myrank .EQ. 0) PRINT*,'# of processors =', np

    !PRINT*,'Myrank =',myrank

    ! get rank coordinates
    CALL MPI_CART_COORDS(comm3d, myrank, ndim, coord, ierr)

    !PRINT*,'Mycoord =',coord


END SUBROUTINE mpi_comm3d_init



SUBROUTINE mpi_get_neighbor_ranks()

    INTEGER :: i, neighbor_coord(3), rank_dest     
    
    ! get neighbor coordinates
        
    ! x- face
    neighbor_coord(:) = (/ coord(1)-1, coord(2), coord(3) /)

    IF( ((neighbor_coord(1) .GT. nranks_x-1) .OR.(neighbor_coord(1) .LT. 0) ) .AND. (.NOT. isperiodic(1))) THEN
        neighbor_rank(1) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(1) = rank_dest
    END IF
       
    ! x+ face
    neighbor_coord(:) = (/ coord(1)+1, coord(2), coord(3) /)

    IF( ((neighbor_coord(1) .GT. nranks_x-1) .OR.(neighbor_coord(1) .LT. 0) ) .AND. (.NOT. isperiodic(1))) THEN
        neighbor_rank(2) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(2) = rank_dest
    END IF
    
    ! y- face
    neighbor_coord(:) = (/ coord(1), coord(2)-1, coord(3) /)

    IF( ((neighbor_coord(2) .GT. nranks_y-1) .OR.(neighbor_coord(2) .LT. 0) ) .AND. (.NOT. isperiodic(2))) THEN
        neighbor_rank(3) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(3) = rank_dest
    END IF
       
    ! y+ face
    neighbor_coord(:) = (/ coord(1), coord(2)+1, coord(3) /)


    IF( ((neighbor_coord(2) .GT. nranks_y-1) .OR.(neighbor_coord(2) .LT. 0) ) .AND. (.NOT. isperiodic(2))) THEN
        neighbor_rank(4) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(4) = rank_dest
    END IF
    
    ! z- face
    neighbor_coord(:) = (/ coord(1), coord(2), coord(3)-1 /)


    IF( ((neighbor_coord(3) .GT. nranks_z-1) .OR.(neighbor_coord(3) .LT. 0) ) .AND. (.NOT. isperiodic(3))) THEN
        neighbor_rank(5) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(5) = rank_dest
    END IF
    
    ! z+ face
    neighbor_coord(:) = (/ coord(1), coord(2), coord(3)+1 /)


    IF( ((neighbor_coord(3) .GT. nranks_z-1) .OR.(neighbor_coord(3) .LT. 0) ) .AND. (.NOT. isperiodic(3))) THEN
        neighbor_rank(6) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(6) = rank_dest
    END IF
    
    
    ! x-edge 1
    neighbor_coord(:) = (/ coord(1), coord(2)-1, coord(3)-1 /)

    IF( ((neighbor_coord(2) .LT. 0) .AND. (.NOT. isperiodic(2))) .AND. &
        ((neighbor_coord(3) .LT. 0) .AND. (.NOT. isperiodic(3)))  ) THEN
        
        neighbor_rank(7) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(7) = rank_dest
    END IF
    
    
    ! x-edge 2
    neighbor_coord(:) = (/ coord(1), coord(2)+1, coord(3)-1 /)

    IF( ((neighbor_coord(2) .GT. nranks_y-1) .AND. (.NOT. isperiodic(2))) .AND. &
        ((neighbor_coord(3) .LT. 0) .AND. (.NOT. isperiodic(3)))  ) THEN
        neighbor_rank(8) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(8) = rank_dest
    END IF
    
    ! x-edge 3
    neighbor_coord(:) = (/ coord(1), coord(2)-1, coord(3)+1 /)

    IF( ((neighbor_coord(2) .LT. 0) .AND. (.NOT. isperiodic(2))) .AND. &
        ((neighbor_coord(3) .GT. nranks_z-1) .AND. (.NOT. isperiodic(3)))  ) THEN
        neighbor_rank(9) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(9) = rank_dest
    END IF
    
    ! x-edge 4
    neighbor_coord(:) = (/ coord(1), coord(2)+1, coord(3)+1 /)
    
    IF( ((neighbor_coord(2) .GT. nranks_y-1) .AND. (.NOT. isperiodic(2))) .AND. &
        ((neighbor_coord(3) .GT. nranks_z-1) .AND. (.NOT. isperiodic(3)))  ) THEN
        neighbor_rank(10) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(10) = rank_dest
    END IF
    
    
     ! y-edge 1
    neighbor_coord(:) = (/ coord(1)-1, coord(2), coord(3)-1 /)

    IF( ((neighbor_coord(1) .LT. 0) .AND. (.NOT. isperiodic(1))) .AND. &
        ((neighbor_coord(3) .LT. 0) .AND. (.NOT. isperiodic(3)))  ) THEN
        neighbor_rank(11) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(11) = rank_dest
    END IF
    
    
    ! y-edge 2
    neighbor_coord(:) = (/ coord(1)+1, coord(2), coord(3)-1 /)

    IF( ((neighbor_coord(1) .GT. nranks_x-1) .AND. (.NOT. isperiodic(1))) .AND. &
        ((neighbor_coord(3) .LT. 0) .AND. (.NOT. isperiodic(3)))  ) THEN
        neighbor_rank(12) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(12) = rank_dest
    END IF
    
    ! y-edge 3
    neighbor_coord(:) = (/ coord(1)-1, coord(2), coord(3)+1 /)

    IF( ((neighbor_coord(1) .LT. 0) .AND. (.NOT. isperiodic(1))) .AND. &
        ((neighbor_coord(3) .GT. nranks_z-1) .AND. (.NOT. isperiodic(3)))  ) THEN
        neighbor_rank(13) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(13) = rank_dest
    END IF
    
    ! y-edge 4
    neighbor_coord(:) = (/ coord(1)+1, coord(2), coord(3)+1 /)

    IF( ((neighbor_coord(1) .GT. nranks_x-1) .AND. (.NOT. isperiodic(1))) .AND. &
        ((neighbor_coord(3) .GT. nranks_z-1) .AND. (.NOT. isperiodic(3)))  ) THEN
        neighbor_rank(14) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(14) = rank_dest
    END IF
    
    
     ! z-edge 1
    neighbor_coord(:) = (/ coord(1)-1, coord(2)-1, coord(3) /)

    IF( ((neighbor_coord(1) .LT. 0) .AND. (.NOT. isperiodic(1))) .AND. &
        ((neighbor_coord(2) .LT. 0) .AND. (.NOT. isperiodic(2)))  ) THEN
        neighbor_rank(15) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(15) = rank_dest
    END IF
    
    
    ! z-edge 2
    neighbor_coord(:) = (/ coord(1)+1, coord(2)-1, coord(3) /)

    IF( ((neighbor_coord(1) .GT. nranks_x-1) .AND. (.NOT. isperiodic(1))) .AND. &
        ((neighbor_coord(2) .LT. 0) .AND. (.NOT. isperiodic(2)))  ) THEN
        neighbor_rank(16) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(16) = rank_dest
    END IF
    
    ! z-edge 3
    neighbor_coord(:) = (/ coord(1)-1, coord(2)+1, coord(3) /)

    IF( ((neighbor_coord(1) .LT. 0) .AND. (.NOT. isperiodic(1))) .AND. &
        ((neighbor_coord(2) .GT. nranks_y-1) .AND. (.NOT. isperiodic(2)))  ) THEN
        neighbor_rank(17) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(17) = rank_dest
    END IF
    
    ! z-edge 4
    neighbor_coord(:) = (/ coord(1)+1, coord(2)+1, coord(3) /)

    IF( ((neighbor_coord(1) .GT. nranks_x-1) .AND. (.NOT. isperiodic(1))) .AND. &
        ((neighbor_coord(2) .GT. nranks_y-1) .AND. (.NOT. isperiodic(2)))  ) THEN
        neighbor_rank(18) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(18) = rank_dest
    END IF
    
    
    ! bottom SW corner
    neighbor_coord(:) = (/ coord(1)-1, coord(2)-1, coord(3)-1 /)

    IF( ((neighbor_coord(1) .LT. 0) .AND. (.NOT. isperiodic(1))) .AND. &
        ((neighbor_coord(2) .LT. 0) .AND. (.NOT. isperiodic(2))) .AND. &
        ((neighbor_coord(3) .LT. 0) .AND. (.NOT. isperiodic(3))) ) THEN
        neighbor_rank(19) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(19) = rank_dest
    END IF
    
    ! bottom SE corner
    neighbor_coord(:) = (/ coord(1)+1, coord(2)-1, coord(3)-1 /)

    IF( ((neighbor_coord(1) .GT. nranks_x-1) .AND. (.NOT. isperiodic(1))) .AND. &
        ((neighbor_coord(2) .LT. 0) .AND. (.NOT. isperiodic(2))) .AND. &
        ((neighbor_coord(3) .LT. 0) .AND. (.NOT. isperiodic(3))) ) THEN
        neighbor_rank(20) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(20) = rank_dest
    END IF
    
    
    ! bottom NW corner
    neighbor_coord(:) = (/ coord(1)-1, coord(2)+1, coord(3)-1 /)

    IF( ((neighbor_coord(1) .LT. 0) .AND. (.NOT. isperiodic(1))) .AND. &
        ((neighbor_coord(2) .GT. nranks_y-1) .AND. (.NOT. isperiodic(2))) .AND. &
        ((neighbor_coord(3) .LT. 0) .AND. (.NOT. isperiodic(3))) ) THEN
        neighbor_rank(21) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(21) = rank_dest
    END IF
    
    
    ! bottom NE corner
        neighbor_coord(:) = (/ coord(1)+1, coord(2)+1, coord(3)-1 /)

    IF( ((neighbor_coord(1) .GT. nranks_x-1) .AND. (.NOT. isperiodic(1))) .AND. &
        ((neighbor_coord(2) .GT. nranks_y-1) .AND. (.NOT. isperiodic(2))) .AND. &
        ((neighbor_coord(3) .LT. 0) .AND. (.NOT. isperiodic(3))) ) THEN
        neighbor_rank(22) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(22) = rank_dest
    END IF
    
    ! top SW corner
    neighbor_coord(:) = (/ coord(1)-1, coord(2)-1, coord(3)+1 /)

    IF( ((neighbor_coord(1) .LT. 0) .AND. (.NOT. isperiodic(1))) .AND. &
        ((neighbor_coord(2) .LT. 0) .AND. (.NOT. isperiodic(2))) .AND. &
        ((neighbor_coord(3) .GT. nranks_z-1) .AND. (.NOT. isperiodic(3))) ) THEN
        neighbor_rank(23) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(23) = rank_dest
    END IF
    
    ! top SE corner
    neighbor_coord(:) = (/ coord(1)+1, coord(2)-1, coord(3)+1 /)

    IF( ((neighbor_coord(1) .GT. nranks_x-1) .AND. (.NOT. isperiodic(1))) .AND. &
        ((neighbor_coord(2) .LT. 0) .AND. (.NOT. isperiodic(2))) .AND. &
        ((neighbor_coord(3) .GT. nranks_z-1) .AND. (.NOT. isperiodic(3))) ) THEN
        neighbor_rank(24) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(24) = rank_dest
    END IF
    
    
    ! top NW corner
    neighbor_coord(:) = (/ coord(1)-1, coord(2)+1, coord(3)+1 /)

    IF( ((neighbor_coord(1) .LT. 0) .AND. (.NOT. isperiodic(1))) .AND. &
        ((neighbor_coord(2) .GT. nranks_y-1) .AND. (.NOT. isperiodic(2))) .AND. &
        ((neighbor_coord(3) .GT. nranks_z-1) .AND. (.NOT. isperiodic(3))) ) THEN
        neighbor_rank(25) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(25) = rank_dest
    END IF
    
    
    ! top NE corner
    neighbor_coord(:) = (/ coord(1)+1, coord(2)+1, coord(3)+1 /)

    IF( ((neighbor_coord(1) .GT. nranks_x-1) .AND. (.NOT. isperiodic(1))) .AND. &
        ((neighbor_coord(2) .GT. nranks_y-1) .AND. (.NOT. isperiodic(2))) .AND. &
        ((neighbor_coord(3) .GT. nranks_z-1) .AND. (.NOT. isperiodic(3))) ) THEN
        neighbor_rank(26) = MPI_PROC_NULL
    ELSE 
        CALL MPI_CART_RANK(comm3d, neighbor_coord, rank_dest, ierr)
        neighbor_rank(26) = rank_dest
    END IF
    
    
    DO i = 0, np-1
    
        CALL MPI_BARRIER(comm3d, ierr)
    
        IF(myrank .EQ. i) THEN
        
            PRINT*,''
            PRINT*,'#######################################################################################'
            WRITE(*,FMT='(" Myrank = ",i2,", Neighbor Faces:  x- = ",i2,", x+ = ",i2,", y- = ",i2,", y+ = ",i2,", z- = ",i2,", z+ = ",i2)') &
            myrank, neighbor_rank(1) ,neighbor_rank(2) ,neighbor_rank(3) ,neighbor_rank(4) ,neighbor_rank(5) ,neighbor_rank(6) 
            WRITE(*,FMT='(" X-edge:  edge1 = ",i2,", edge2 = ",i2,", edge3 = ",i2,", edge4 = ",i2)') &
            neighbor_rank(7) ,neighbor_rank(8) ,neighbor_rank(9) ,neighbor_rank(10) 
            WRITE(*,FMT='(" Y-edge:  edge1 = ",i2,", edge2 = ",i2,", edge3 = ",i2,", edge4 = ",i2)') &
            neighbor_rank(11) ,neighbor_rank(12) ,neighbor_rank(13) ,neighbor_rank(13) 
            WRITE(*,FMT='(" Z-edge:  edge1 = ",i2,", edge2 = ",i2,", edge3 = ",i2,", edge4 = ",i2)') &
            neighbor_rank(15) ,neighbor_rank(16) ,neighbor_rank(17) ,neighbor_rank(18)
            WRITE(*,FMT='(" Bottom-corners:  SW = ",i2,", SE = ",i2,", NW = ",i2,", NE = ",i2)') &
            neighbor_rank(19) ,neighbor_rank(20) ,neighbor_rank(21) ,neighbor_rank(22)
            WRITE(*,FMT='(" Top-corners:  SW = ",i2,", SE = ",i2,", NW = ",i2,", NE = ",i2)') &
            neighbor_rank(23) ,neighbor_rank(24) ,neighbor_rank(25) ,neighbor_rank(26)
            
            PRINT*,'My coordinates =',coord
            PRINT*,'######################################################################################'
            PRINT*,''
        
        END IF

        CALL MPI_BARRIER(comm3d, ierr)
    
    END DO
    
    
    
    
END SUBROUTINE mpi_get_neighbor_ranks



SUBROUTINE rma_init()

    INTEGER :: i, ierr, num
    
  
    ! determine array element size in bytes
    CALL MPI_SIZEOF(DBLE(1), sizedouble, ierr)
    num = 2 * 26 + 2 * (nx*ny + ny*nz + nz*nx) * nb + 4 * (nx + ny + nz) * nb * nb + 8 * nb * nb
    winsize =  num * sizedouble
    
    IF(myrank .EQ. 0) PRINT*,'Winsize (Mbytes) = ',winsize*1e-6 
    
    ! create local rma memory window and allocate memory block
    CALL MPI_WIN_ALLOCATE(winsize, sizedouble, MPI_INFO_NULL, comm3d, rma_cmem, win, ierr)
    CALL C_F_POINTER(rma_cmem, A, (/ num /))


    ! Note: The RMA window is broken up into 26 chunks corresponding to ghost cells at the 
    ! * 6 cell face boundaries: x-, x+, y-, y+, z-, z+ (in that order)
    ! * 12 cell edge  boundaries, i.e. 4 edges along to each coordinate direction: x(1,2,3,4), y(1,2,3,4), z(1,2,3,4) 
    ! * 8 cell corners, bottom (SW, SE, NW, NE), top(SW, SE, NW, NE)
    ! The first two slots in each chunk are reserved for the data recv/sent signals from the neighbor corresponding to that boundary.
    
    A(:) = 0.d0

    ! set rma window displacement units for the different chunks
    disp(1)  = 0
    disp(2)  = disp(1)  + nb * (ny*nz) + 2 
    disp(3)  = disp(2)  + nb * (ny*nz) + 2
    disp(4)  = disp(3)  + nb * (nx*nz) + 2
    disp(5)  = disp(4)  + nb * (nx*nz) + 2
    disp(6)  = disp(5)  + nb * (nx*ny) + 2
    disp(7)  = disp(6)  + nb * (nx*ny) + 2
    disp(8)  = disp(7)  + nx * nb * nb + 2
    disp(9)  = disp(8)  + nx * nb * nb + 2
    disp(10) = disp(9)  + nx * nb * nb + 2
    disp(11) = disp(10) + nx * nb * nb + 2
    disp(12) = disp(11) + ny * nb * nb + 2
    disp(13) = disp(12) + ny * nb * nb + 2
    disp(14) = disp(13) + ny * nb * nb + 2
    disp(15) = disp(14) + ny * nb * nb + 2
    disp(16) = disp(15) + nz * nb * nb + 2
    disp(17) = disp(16) + nz * nb * nb + 2
    disp(18) = disp(17) + nz * nb * nb + 2
    disp(19) = disp(18) + nz * nb * nb + 2
    disp(20) = disp(19) + nb * nb + 2
    disp(21) = disp(20) + nb * nb + 2
    disp(22) = disp(21) + nb * nb + 2
    disp(23) = disp(22) + nb * nb + 2
    disp(24) = disp(23) + nb * nb + 2
    disp(25) = disp(24) + nb * nb + 2
    disp(26) = disp(25) + nb * nb + 2
       
    ! shuffled displacement units (for convenience later when doing PUTs into neighbor boundaries)  
    shuffled_disp = (/ disp(2), disp(1), disp(4), disp(3), disp(6), disp(5),  &
                       disp(10), disp(9), disp(8), disp(7),     &
                       disp(14), disp(13), disp(12), disp(11),  & 
                       disp(18), disp(17), disp(16), disp(15),  &
                       disp(26), disp(25), disp(24), disp(23),  &
                       disp(22), disp(21), disp(20), disp(19)     /)


    ! sizes of PUT buffers corresponding to the 26 boundaries 
    send_buffer_size = (/ ny*nz*nb, ny*nz*nb, nx*nz*nb, nx*nz*nb, nx*ny*nb, nx*ny*nb, &
                           nx*nb*nb, nx*nb*nb, nx*nb*nb, nx*nb*nb, &
                           ny*nb*nb, ny*nb*nb, ny*nb*nb, ny*nb*nb, &
                           nz*nb*nb, nz*nb*nb, nz*nb*nb, nz*nb*nb, &
                           nb*nb, nb*nb, nb*nb, nb*nb, nb*nb, nb*nb, nb*nb, nb*nb /)

    !PRINT*,'MPI initialization completed.'

END SUBROUTINE rma_init


! shut down MPI
SUBROUTINE close_rma_window()

    INTEGER :: ierr

    ! end rma exposure epoch
    CALL MPI_WIN_UNLOCK(neighbor, win, ierr)
    
    ! free rma window
    CALL MPI_WIN_FREE(win, ierr)
    

    NULLIFY(A)

END SUBROUTINE close_rma_window


! start MPI passive RMA epoch
SUBROUTINE open_exposure_epoch()

    INTEGER :: ierr
    
    !PRINT*,'Opening passive exposure epoch.'
    
    ! start exposure epoch so that neighbors can access my rma window
    CALL MPI_WIN_LOCK(MPI_LOCK_SHARED, neighbor, MPI_MODE_NOCHECK, win, ierr)
    
    !PRINT*,'Passive exposure epoch now open.'
    
END SUBROUTINE open_exposure_epoch


SUBROUTINE print_array(offset)

    INTEGER, INTENT(IN) :: offset
    INTEGER :: n

    DO n = 0, np-1

    CALL MPI_BARRIER(comm3d, ierr)

        IF(myrank .EQ. n) THEN 
    
    
            PRINT*,''
            PRINT*,'Myrank = ',myrank
            PRINT*,''
            PRINT*,'u(x- face) ='
            DO k = nz+nb, 1-nb, -1
            DO j = 1-nb, ny+nb
                WRITE(*,FMT='(f5.1)', ADVANCE='NO') u(1-offset,j,k)
            END DO    
            PRINT*,''
            END DO
    
            PRINT*,''
            PRINT*,'u(x+ face) = '
            DO k = nz+nb, 1-nb, -1
            DO j = 1-nb, ny+nb
                WRITE(*,FMT='(f5.1)', ADVANCE='NO') u(nx+offset,j,k)
            END DO    
            PRINT*,''
            END DO

            !PRINT*,''
            !PRINT*,'u(y- face) ='
            !DO k = nz+nb, 1-nb, -1
            !DO i = 1-nb, nx+nb
            !    WRITE(*,FMT='(f5.1)', ADVANCE='NO') u(i,1-nb,k)
            !END DO    
            !PRINT*,''
            !END DO
    
            !PRINT*,''
            !PRINT*,'u(y+ face) ='
            !  DO k = nz+nb, 1-nb, -1
            !DO i = 1-nb, nx+nb
            !    WRITE(*,FMT='(f5.1)', ADVANCE='NO') u(i,ny+nb,k)
            !END DO    
            !PRINT*,''
            !END DO
    
            PRINT*,''
            PRINT*,'u(z- face) ='
            DO j = ny+nb, 1-nb, -1
            DO i = 1-nb, nx+nb
                WRITE(*,FMT='(f5.1)', ADVANCE='NO') u(i,j,1-offset)
            END DO    
            PRINT*,''
            END DO
    
            PRINT*,''
            PRINT*,'u(z+ face) ='
            DO j = ny+nb, 1-nb, -1
            DO i = 1-nb, nx+nb
                WRITE(*,FMT='(f5.1)', ADVANCE='NO') u(i,j,nz+offset)
            END DO    
            PRINT*,''
            END DO
    
            CALL sleep(1)
        
        END IF

        CALL MPI_BARRIER(comm3d, ierr)
    
    END DO


END SUBROUTINE print_array



END PROGRAM test7_onesided