MODULE constants_mod

IMPLICIT NONE


!********************************
! simulation and grid parameters
!********************************
INTEGER, PARAMETER :: nx       = 64    ! grid size per MPI rank
INTEGER, PARAMETER :: ny       = 64    	     
INTEGER, PARAMETER :: nz       = 64 	     
INTEGER, PARAMETER :: nb       = 6    ! number of boundary cells (NEEDS TO BE >= 5)

REAL*8, PARAMETER  :: t_end = 2.5d0         ! simulation end time
INTEGER, PARAMETER :: maxsteps = 1000000000   ! max. number of time steps
REAL*8, PARAMETER  :: COUR_def = 0.35d0       ! Default Courant number (needs to be <= 0.5, prerferably <=0.35 if turbulence driver is on)

LOGICAL, PARAMETER :: print_debug     = .FALSE.
LOGICAL, PARAMETER :: override_checks = .TRUE.

REAL*8,  PARAMETER :: FOURPI = 16.d0*ATAN(1.d0)
REAL*8,  PARAMETER :: TWOPI  = 8.d0*ATAN(1.d0)
INTEGER, PARAMETER :: s1_def = 1234, s2_def = 5678, s3_def = 9123  ! default RNG seed values

!*************************************
! MPI domain decomposition parameters
!*************************************
INTEGER, PARAMETER :: nranks_x = 1  ! # of ranks along x direction
INTEGER, PARAMETER :: nranks_y = 1  ! # of ranks along y direction
INTEGER, PARAMETER :: nranks_z = 1  ! # of ranks along z direction

!********************
! physics parameters
!********************
REAL*8 :: sound_speed = 1.d0  ! constant isothermal sound speed
INTEGER, PARAMETER :: boundary_type = 2  ! 1: open , 2: periodic

!*************************
!Turbulence parameters
!*************************
LOGICAL, PARAMETER :: drive_turbulence  = .FALSE.  ! turbulence driver switch
REAL*8, PARAMETER :: dt_turb = 0.05d0
INTEGER, PARAMETER :: k_max = 2    ! Largest wave number for the force. NEEDS TO BE < MIN(nx,ny,nz)/2
REAL*8, PARAMETER :: Edot = 1.d-3  ! Specify the constant rate of kinetic energy injection    
REAL*8, PARAMETER :: Bx_0 = 1.d-5   ! Uniform initial magnetic field
REAL*8, PARAMETER :: By_0 = 0.d0
REAL*8, PARAMETER :: Bz_0 = 0.d0

!***************
! MPI variables
!***************
INTEGER :: myrank, my_coord(3)

!**********************
! MHD Solver parameters
!**********************
INTEGER, PARAMETER :: riemann_solver_type = 2 ! 1: ROE, 2: HLLD, 3: HLLE
INTEGER, PARAMETER :: mhd_solver_type = 2 ! 1 = TVD, 2 = Van Leer

!***********************************
! Passive variable solver parameters
!***********************************
INTEGER, PARAMETER :: npass = 1

!********************
! Restart Parameters
!********************
LOGICAL, PARAMETER :: read_from_restart_file = .FALSE.
REAL*8, PARAMETER :: restart_frequency = 0.5d0

!************************
! File output Parameters
!************************
REAL*8, PARAMETER :: dump_frequency = 0.02d0
CHARACTER(LEN=300), PARAMETER :: output_filepath = 'Output'     
        
LOGICAL, PARAMETER :: parallel_io = .FALSE.
INTEGER, PARAMETER :: parallel_filesize = nx * ny * nz * (7 + npass)  

LOGICAL, PARAMETER :: double_prec = .FALSE.

LOGICAL, PARAMETER :: output_xcut = .FALSE.
LOGICAL, PARAMETER :: output_ycut = .FALSE.
LOGICAL, PARAMETER :: output_zcut = .FALSE.
LOGICAL, PARAMETER :: output_xydiagcut = .FALSE.
LOGICAL, PARAMETER :: output_yzdiagcut = .FALSE.
LOGICAL, PARAMETER :: output_xy_plane = .FALSE.
LOGICAL, PARAMETER :: output_xyz_cube = .TRUE.

END MODULE constants_mod