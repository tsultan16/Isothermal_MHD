MODULE constants_mod

IMPLICIT NONE


INTEGER, PARAMETER :: nranks_test = 1

!********************************
! simulation and grid parameters
!********************************
INTEGER, PARAMETER :: nx       = 30    ! grid size per MPI rank
INTEGER, PARAMETER :: ny       = 30    	     
INTEGER, PARAMETER :: nz       = 30   	     
INTEGER, PARAMETER :: nb       = 5     ! number of boundary cells (needs to be >= 5)

INTEGER, PARAMETER :: maxsteps = 5    ! max. number of time steps
INTEGER, PARAMETER :: tskip    = 2     ! file output frequency
REAL*8, PARAMETER  :: COUR = 0.5D0     ! Courant number (needs to be <= 0.5)

LOGICAL, PARAMETER :: print_debug       = .FALSE.
LOGICAL, PARAMETER :: override_checks   = .TRUE.

REAL*8,  PARAMETER :: FOURPI  = 16.d0*ATAN(1.d0)
REAL*8 , PARAMETER :: third  = 1.d0 / 3.d0 

!*************************************
! MPI domain decomposition parameters
!*************************************
INTEGER, PARAMETER :: nranks_x = 3  ! # of ranks along x direction
INTEGER, PARAMETER :: nranks_y = 3  ! # of ranks along y direction
INTEGER, PARAMETER :: nranks_z = 3  ! # of ranks along z direction

!********************
! physics parameters
!********************
REAL*8 :: sound_speed = 1.d0  ! constant isothermal sound speed
INTEGER, PARAMETER :: boundary_type = 2  ! 1: open , 2: periodic

!***************
! MPI variables
!***************
INTEGER :: myrank, my_coord(3)

!*****************
! init parameters
!*****************

!********************
! Restart Parameters
!********************
LOGICAL, PARAMETER :: read_from_restart_file = .FALSE.


!************************
! File output Parameters
!************************
CHARACTER(LEN=300), PARAMETER :: output_filepath = 'Output/Snapshots/'  ! '/export/data/local/tanzid/Output/Snapshots/'  
        
LOGICAL, PARAMETER :: parallel_io = .TRUE.
INTEGER, PARAMETER :: parallel_filesize = nx * ny * nz * 10 

LOGICAL, PARAMETER :: output_xcut = .FALSE.
LOGICAL, PARAMETER :: output_ycut = .FALSE.
LOGICAL, PARAMETER :: output_zcut = .FALSE.
LOGICAL, PARAMETER :: output_xydiagcut = .FALSE.
LOGICAL, PARAMETER :: output_yzdiagcut = .FALSE.
LOGICAL, PARAMETER :: output_xy_plane = .FALSE.
LOGICAL, PARAMETER :: output_xyz_cube = .TRUE.

END MODULE constants_mod