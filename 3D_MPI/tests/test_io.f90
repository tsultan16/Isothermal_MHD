!!!!
PROGRAM test_io


IMPLICIT NONE

INCLUDE 'mpif.h' 


REAL*8, ALLOCATABLE :: array(:)
INTEGER :: byte_size, arr_size
REAL*8 :: t1, t2

! set array size in bytes
byte_size = 2000000
arr_size = byte_size/8

ALLOCATE(array(arr_size))

array = 1.0

t1 = MPI_Wtime()

OPEN(UNIT=10, FILE = 'output.dat', FORM = 'UNFORMATTED', ACCESS = 'SEQUENTIAL')


WRITE(10) array


CLOSE(UNIT=10)

t2 = MPI_Wtime()


PRINT*,''
PRINT*,'File size (Mb) =',byte_size*1e-6
PRINT*,'File I/O time (sec) = ',t2-t1
PRINT*,''

DEALLOCATE(array)

END PROGRAM test_io