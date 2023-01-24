MODULE constants_mod

IMPLICIT NONE


!********************************
! simulation and grid parameters
!********************************
INTEGER, PARAMETER :: nx       = 64    ! grid size per MPI rank
INTEGER, PARAMETER :: ny       = 64    	     
INTEGER, PARAMETER :: nz       = 64  	     
INTEGER, PARAMETER :: nb       = 5     ! number of boundary cells (NEEDS TO BE >= 5)

REAL*8, PARAMETER  :: t_end = 10.d0           ! simulation end time


LOGICAL, PARAMETER :: print_debug     = .FALSE.
LOGICAL, PARAMETER :: override_checks = .TRUE.

REAL*4,  PARAMETER :: FOURPI = 16.d0*ATAN(1.d0)
REAL*4,  PARAMETER :: TWOPI  = 8.d0*ATAN(1.d0)


!*************************************
! MPI domain decomposition parameters
!*************************************
INTEGER, PARAMETER :: nranks_x = 4  ! # of ranks along x direction
INTEGER, PARAMETER :: nranks_y = 4  ! # of ranks along y direction
INTEGER, PARAMETER :: nranks_z = 4  ! # of ranks along z direction

!********************
! physics parameters
!********************
REAL*4 :: sound_speed = 1.d0  ! constant isothermal sound speed

!***********************************
! Passive variable solver parameters
!***********************************
INTEGER, PARAMETER :: npass = 1

!************************
! File output Parameters
!************************
CHARACTER(LEN=300), PARAMETER :: output_filepath = '/data/uchu/tanzid/Isothermal_MHD/3D_MPI/Output'  
REAL*4, PARAMETER :: dump_frequency = 0.2         


END MODULE constants_mod