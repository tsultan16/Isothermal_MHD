MODULE constants_mod

IMPLICIT NONE


! simulation and grid parameters
INTEGER, PARAMETER :: nx       = 16      	     ! mesh size
INTEGER, PARAMETER :: ny       = 64    	     
INTEGER, PARAMETER :: nz       = 64    	     
INTEGER, PARAMETER :: nb       = 5               ! number of boundary cells

INTEGER, PARAMETER :: maxsteps = 100             ! max. number of time steps
INTEGER, PARAMETER :: tskip    = 5             ! file output frequency

LOGICAL, PARAMETER :: print_debug       = .FALSE.
LOGICAL, PARAMETER :: override_checks   = .TRUE.

REAL*8,  PARAMETER :: FOURPI  = 16.d0*ATAN(1.d0)
REAL*8 , PARAMETER :: third  = 1.d0 / 3.d0 


! physics parameters
REAL*8 :: sound_speed = 1.d0  ! constant isothermal sound speed


! init parameters


! File output parameters



!MPI Parameters
INTEGER, PARAMETER :: ndims  = 1 ! domain decomposition dimensions
INTEGER, PARAMETER :: nvars_particles = 7 ! number of variables per particle (x,y,z,ux,uy,uz)
INTEGER, PARAMETER :: nvars_fields = 3  ! number of field variables (E/B/J)
INTEGER, PARAMETER :: nranks_x = 2
INTEGER, PARAMETER :: nranks_y = 1
INTEGER, PARAMETER :: nranks_z = 1 



END MODULE constants_mod