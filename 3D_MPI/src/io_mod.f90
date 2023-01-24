MODULE io_mod

USE constants_mod
USE grid_data_mod
USE mpi_domain_mod

IMPLICIT NONE

CONTAINS


! This subroutine writes the MHD state variables grid array from a single time-step 
! into an unformatted (i.e. binary) file, i.e. a stream of bytes.
! Each MPI rank writes to it's on file.
SUBROUTINE writetofile_unformatted(t_dump)


    INTEGER, INTENT(IN) :: t_dump
	
    INTEGER :: i, j, k
    REAL*8 :: dx, vx, vy, vz, vpar, vperp, bpar, bperp, isqr2
    REAL*8 :: x, y, z
    CHARACTER(LEN=300) :: filename1, filename2, filename3, filename4, filename5, filename6, filename7
    CHARACTER(LEN=6) :: uniti, coordx, coordy, coordz
    CHARACTER(LEN=20) :: riemann_solver_prefix
    
    
    IF(t_dump<10) THEN
        WRITE(uniti,'(I1.1)') t_dump
    ELSE IF(t_dump>=10 .and. t_dump<100) THEN
        WRITE(uniti,'(I2.2)') t_dump
    ELSE IF(t_dump>=100 .and. t_dump<1000) THEN
        WRITE (uniti,'(I3.3)') t_dump
    ELSE IF(t_dump>=1000 .and. t_dump<10000) THEN
        WRITE (uniti,'(I4.3)') t_dump
    ELSE IF(t_dump>=10000 .and. t_dump<100000) THEN
        WRITE (uniti,'(I5.3)') t_dump  
    END IF
      
    
    WRITE(coordx,'(I1.1)') my_coord(1)
    WRITE(coordy,'(I1.1)') my_coord(2)
    WRITE(coordz,'(I1.1)') my_coord(3)
    
    
    IF(riemann_solver_type .EQ. 1) THEN
        riemann_solver_prefix = 'ROE'
    ELSE IF(riemann_solver_type .EQ. 2) THEN
        riemann_solver_prefix = 'HLLD'
    ELSE IF(riemann_solver_type .EQ. 3) THEN
        riemann_solver_prefix = 'HLLE'
    END IF
        
    
    filename1 = TRIM(output_filepath)//TRIM('/Snapshots/')//TRIM('xcut_')//TRIM(riemann_solver_prefix)//('_rank=')//TRIM(coordx)//TRIM(coordy) &
                //TRIM(coordz)//TRIM('_dump=')//TRIM(uniti)//TRIM('.dat')
                
    filename2 = TRIM(output_filepath)//TRIM('/Snapshots/')//TRIM('ycut_rank=')//TRIM(coordx)//TRIM(coordy) &
                //TRIM(coordz)//TRIM('_dump=')//TRIM(uniti)//TRIM('.dat')
                
    filename3 = TRIM(output_filepath)//TRIM('/Snapshots/')//TRIM('zcut_rank=')//TRIM(coordx)//TRIM(coordy) & 
                //TRIM(coordz)//TRIM('_dump=')//TRIM(uniti)//TRIM('.dat')
                
    filename4 = TRIM(output_filepath)//TRIM('/Snapshots/')//TRIM('xydiagcut_rank=')//TRIM(coordx)//TRIM(coordy) &
                //TRIM(coordz)//TRIM('_dump=')//TRIM(uniti)//TRIM('.dat')
                
    filename5 = TRIM(output_filepath)//TRIM('/Snapshots/')//TRIM('yzdiagcut_rank=')//TRIM(coordx)//TRIM(coordy) & 
                //TRIM(coordz)//TRIM('_dump=')//TRIM(uniti)//TRIM('.dat')

    filename6 = TRIM(output_filepath)//TRIM('/Snapshots/')//TRIM('xyplane_rank=')//TRIM(coordx)//TRIM(coordy) & 
                //TRIM(coordz)//TRIM('_dump=')//TRIM(uniti)//TRIM('.dat')
 
    filename7 = TRIM(output_filepath)//TRIM('/Snapshots/')//TRIM('xyzcube_rank=')//TRIM(coordx)//TRIM(coordy) & 
                //TRIM(coordz)//TRIM('_dump=')//TRIM(uniti)//TRIM('.dat')
                
    dx = 1.d0 / DBLE(MAX(nx*nranks_x,ny*nranks_y,nz*nranks_z))
    isqr2 = 1.d0 / SQRT(2.d0)
 

    IF(output_xcut) THEN
 
    OPEN(UNIT=11,FILE=filename1, FORM = 'UNFORMATTED', ACCESS = 'STREAM')
    
    
    j = ny/2!ny-2
    k = nz/2
    
    DO i =1,nx
        vx = q_3d(i,j,k,2)/q_3d(i,j,k,1)
        vy = q_3d(i,j,k,3)/q_3d(i,j,k,1)
        vz = q_3d(i,j,k,4)/q_3d(i,j,k,1)
        
        IF(double_prec) THEN
            WRITE(11) q_3d(i,j,k,1),vx,vy,vz,q_3d(i,j,k,5),q_3d(i,j,k,6),q_3d(i,j,k,7)
        ELSE
            WRITE(11) SNGL(q_3d(i,j,k,1)),SNGL(vx),SNGL(vy),SNGL(vz),SNGL(q_3d(i,j,k,5)),SNGL(q_3d(i,j,k,6)),SNGL(q_3d(i,j,k,7))
        END IF        
    END DO

    CLOSE(UNIT=11)

    END IF


    IF(output_ycut) THEN

    OPEN(UNIT=12,FILE=filename2, FORM = 'UNFORMATTED', ACCESS = 'STREAM')

    i = nx/2
    k = nz/2
    
    DO j =1,ny
        vx = q_3d(i,j,k,2)/q_3d(i,j,k,1)
        vy = q_3d(i,j,k,3)/q_3d(i,j,k,1)
        vz = q_3d(i,j,k,4)/q_3d(i,j,k,1)
        
        IF(double_prec) THEN
            WRITE(12) q_3d(i,j,k,1),vx,vy,vz,q_3d(i,j,k,5),q_3d(i,j,k,6),q_3d(i,j,k,7)
        ELSE
            WRITE(12) SNGL(q_3d(i,j,k,1)),SNGL(vx),SNGL(vy),SNGL(vz),SNGL(q_3d(i,j,k,5)),SNGL(q_3d(i,j,k,6)),SNGL(q_3d(i,j,k,7))
        END IF  
    END DO

    CLOSE(UNIT=12)

    END IF
    
    IF(output_zcut) THEN

    OPEN(UNIT=13,FILE=filename3, FORM = 'UNFORMATTED', ACCESS = 'STREAM')

    i = nx/2
    j = ny/2
    
    DO k =1,nz
        vx = q_3d(i,j,k,2)/q_3d(i,j,k,1)
        vy = q_3d(i,j,k,3)/q_3d(i,j,k,1)
        vz = q_3d(i,j,k,4)/q_3d(i,j,k,1)
        
        IF(double_prec) THEN
            WRITE(13) q_3d(i,j,k,1),vx,vy,vz,q_3d(i,j,k,5),q_3d(i,j,k,6),q_3d(i,j,k,7)
        ELSE
            WRITE(13) SNGL(q_3d(i,j,k,1)),SNGL(vx),SNGL(vy),SNGL(vz),SNGL(q_3d(i,j,k,5)),SNGL(q_3d(i,j,k,6)),SNGL(q_3d(i,j,k,7))
        END IF   
        
    END DO

    CLOSE(UNIT=13)

    END IF


    IF(output_xydiagcut) THEN

    OPEN(UNIT=14,FILE=filename4, FORM = 'UNFORMATTED', ACCESS = 'STREAM')

    k = nz/2
    
    DO j =1,ny
    DO i =1,nx
        
        IF(i .EQ. j) THEN 

            vx = q_3d(i,j,k,2)/q_3d(i,j,k,1)
            vy = q_3d(i,j,k,3)/q_3d(i,j,k,1)
            vz = q_3d(i,j,k,4)/q_3d(i,j,k,1)

            vperp = isqr2 * (vy + vx) 
            vpar  = isqr2 * (vy - vx)
            bperp = isqr2 * (q_3d(i,j,k,6) + q_3d(i,j,k,5))
            bpar  = isqr2 * (q_3d(i,j,k,6) - q_3d(i,j,k,5))            
        
            IF(double_prec) THEN
                WRITE(14) q_3d(i,j,k,1),vperp,vpar,vz,bperp,bpar,q_3d(i,j,k,7)
            ELSE
                WRITE(14) SNGL(q_3d(i,j,k,1)),SNGL(vperp),SNGL(vpar),SNGL(vz),SNGL(bperp),SNGL(bpar),SNGL(q_3d(i,j,k,7))
            END IF  
        END IF      
        
    END DO
    END DO

    CLOSE(UNIT=14)

    END IF


    IF(output_yzdiagcut) THEN

    OPEN(UNIT=15,FILE=filename5, FORM = 'UNFORMATTED', ACCESS = 'STREAM')
    
    i = nx/2
    
    DO k =1,nz
    DO j =1,ny
        
        IF(j .EQ. k) THEN 

            vx = q_3d(i,j,k,2)/q_3d(i,j,k,1)
            vy = q_3d(i,j,k,3)/q_3d(i,j,k,1)
            vz = q_3d(i,j,k,4)/q_3d(i,j,k,1)

            vperp = isqr2 * (vz + vy) 
            vpar  = isqr2 * (vz - vy)
            bperp = isqr2 * (q_3d(i,j,k,7) + q_3d(i,j,k,6))
            bpar  = isqr2 * (q_3d(i,j,k,7) - q_3d(i,j,k,6))            
        
            IF(double_prec) THEN
                WRITE(15) q_3d(i,j,k,1),vperp,vpar,vz,bperp,bpar,q_3d(i,j,k,7)
            ELSE
                WRITE(15) SNGL(q_3d(i,j,k,1)),SNGL(vperp),SNGL(vpar),SNGL(vz),SNGL(bperp),SNGL(bpar),SNGL(q_3d(i,j,k,7))
            END IF  
        END IF      
        
    END DO
    END DO
    
    CLOSE(UNIT=15)
    
    END IF

    IF(output_xy_plane) THEN

    OPEN(UNIT=16,FILE=filename6, FORM = 'UNFORMATTED', ACCESS = 'STREAM')
    
    k = nz/2
    
    DO j =1,ny
    DO i =1,nx
        
        vx = q_3d(i,j,k,2)/q_3d(i,j,k,1)
        vy = q_3d(i,j,k,3)/q_3d(i,j,k,1)
        vz = q_3d(i,j,k,4)/q_3d(i,j,k,1)
        
        IF(double_prec) THEN
            WRITE(16) q_3d(i,j,k,1),vx,vy,vz,q_3d(i,j,k,5),q_3d(i,j,k,6),q_3d(i,j,k,7)
        ELSE
            WRITE(16) SNGL(q_3d(i,j,k,1)),SNGL(vx),SNGL(vy),SNGL(vz),SNGL(q_3d(i,j,k,5)),SNGL(q_3d(i,j,k,6)),SNGL(q_3d(i,j,k,7)) 
       END IF  
       
    END DO
    END DO
    
    CLOSE(UNIT=16)
    
    END IF
    

    IF(output_xyz_cube) THEN

    OPEN(UNIT=17,FILE=filename7, FORM = 'UNFORMATTED', ACCESS = 'STREAM')
    
    DO k =1,nz 
    DO j =1,ny
    DO i =1,nx
        IF(npass .GT. 0) THEN
            IF(double_prec) THEN
                WRITE(17) q_3d(i,j,k,1),q_3d(i,j,k,2),q_3d(i,j,k,3),q_3d(i,j,k,4),q_3d(i,j,k,5),q_3d(i,j,k,6),q_3d(i,j,k,7),qpass_3d(i,j,k,1)
            ELSE
                WRITE(17) SNGL(q_3d(i,j,k,1)),SNGL(q_3d(i,j,k,2)),SNGL(q_3d(i,j,k,3)),SNGL(q_3d(i,j,k,4)),SNGL(q_3d(i,j,k,5)),SNGL(q_3d(i,j,k,6)), &
                          SNGL(q_3d(i,j,k,7)),SNGL(qpass_3d(i,j,k,1))
            END IF  
       ELSE
            IF(double_prec) THEN
                WRITE(17) q_3d(i,j,k,1),q_3d(i,j,k,2),q_3d(i,j,k,3),q_3d(i,j,k,4),q_3d(i,j,k,5),q_3d(i,j,k,6),q_3d(i,j,k,7)
            ELSE
                WRITE(17) SNGL(q_3d(i,j,k,1)),SNGL(q_3d(i,j,k,2)),SNGL(q_3d(i,j,k,3)),SNGL(q_3d(i,j,k,4)),SNGL(q_3d(i,j,k,5)),SNGL(q_3d(i,j,k,6)), &
                          SNGL(q_3d(i,j,k,7))
            END IF
       END IF
    END DO    
    END DO
    END DO
    
    CLOSE(UNIT=17)
    
    END IF
    
    
END SUBROUTINE writetofile_unformatted


! Parallel version of file output routine. Each mpi rank writes
! to a different portion of a common file visible to all ranks.
! Use this only when a high-performance parallel filesystem is 
! available (e.g. Lustre, GPFS, etc.), otherwise this might be slower
! than each rank writing to it's own file (probably not a good idea 
! to use this routine if running on a single node).
!
! The file is partitioned into (equally-sized) contiguous blocks:
!
!        |---P_O---|---P_1---|  ....  |---P_N-1---|
!        ^                                        ^
!  head of file                                end of file
!
! where |---P_ix---| is the block for the "ix_th" process (ix = 0,1,...,n_processes-1)
!
! In our 3d domain decomposition, each rank has an x,y,z coordinate. We flattened 
! this coordinate along the x-direction, so that the flattening index is 
! ix = my_coord(1) + nranks_x * my_coord(2) + nranks_x * nranks_y * my_coord(3),
! and blocks are assigned in the order of increasing flattened index.
SUBROUTINE writetofile_unformatted_parallel(t_dump, parallel_io_req, file_handle, initial, final, write_status)


    INTEGER, INTENT(IN) :: t_dump
    INTEGER, INTENT(INOUT) :: parallel_io_req, file_handle
	LOGICAL, INTENT(IN) :: initial, final
    LOGICAL, INTENT(OUT) :: write_status
    INTEGER :: ierr, status(MPI_STATUS_SIZE), bytesize
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset
    CHARACTER(LEN=300) :: filename
    CHARACTER(LEN=50) :: uniti

    IF(double_prec) THEN
        bytesize = 8
    ELSE
        bytesize = 4
    END IF
            
    IF(t_dump<10) THEN
        WRITE(uniti,'(I1.1)') t_dump
    ELSE IF(t_dump>=10 .and. t_dump<100) THEN
        WRITE(uniti,'(I2.2)') t_dump
    ELSE IF(t_dump>=100 .and. t_dump<1000) THEN
        WRITE (uniti,'(I3.3)') t_dump
    ELSE IF(t_dump>=1000 .and. t_dump<10000) THEN
        WRITE (uniti,'(I4.3)') t_dump
    ELSE IF(t_dump>=10000 .and. t_dump<100000) THEN
        WRITE (uniti,'(I5.3)') t_dump  
    END IF
      
    ! wait for the previous file write to complete, then close the old file  
    IF(.NOT. initial) THEN
        CALL MPI_WAIT(parallel_io_req, status, ierr)
        CALL MPI_FILE_CLOSE(file_handle, ierr)
    END IF    
      
    ! now start the next non-blocking file write
    IF(.NOT. final) THEN  
        
        filename = TRIM(output_filepath)//TRIM('/Snapshots/')//TRIM('mhd_state_parallel_dump=')//TRIM(uniti)//TRIM('.dat')    
        
        ! set the file offset (i.e. number of bytes to skip from the start of the file) 
        offset = ( my_coord(1) + nranks_x * my_coord(2) + nranks_x * nranks_y * my_coord(3) ) * parallel_filesize * bytesize  

        ! open new parallel file
        CALL MPI_FILE_OPEN(comm3d, filename, MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, file_handle, ierr)

        ! write to parallel file (Note: this is a non-blocking MPI operation)
        IF(double_prec) THEN
            CALL MPI_FILE_IWRITE_AT(file_handle, offset, file_buffer_dble, parallel_filesize, MPI_DOUBLE_PRECISION, parallel_io_req, ierr)
        ELSE
            CALL MPI_FILE_IWRITE_AT(file_handle, offset, file_buffer_sngl, parallel_filesize, MPI_REAL, parallel_io_req, ierr)
        END IF
        
        IF(ierr .NE. MPI_SUCCESS) THEN
            write_status = .FALSE.
        ELSE 
            write_status = .TRUE.        
        END IF

    END IF

END SUBROUTINE writetofile_unformatted_parallel


! this subroutine copies mhd state variables into our parallel file io buffer
SUBROUTINE load_parallel_file_buffer(t_sim)
	
    REAL*8, INTENT(IN) :: t_sim
    INTEGER :: i, j, k, ix, ipass
    REAL*8 :: x, y, z, dx
 
    dx = 1.d0 / DBLE(MAX(nx*nranks_x,ny*nranks_y,nz*nranks_z))
  
    ix = 1

    IF(double_prec) THEN

        ! load up output data into the io buffer 
        DO k = 1, nz    
        DO j = 1, ny
        DO i = 1, nx
            
            file_buffer_dble(ix+0) = q_3d(i,j,k,1)
            file_buffer_dble(ix+1) = q_3d(i,j,k,2)
            file_buffer_dble(ix+2) = q_3d(i,j,k,3)
            file_buffer_dble(ix+3) = q_3d(i,j,k,4)
            file_buffer_dble(ix+4) = q_3d(i,j,k,5)
            file_buffer_dble(ix+5) = q_3d(i,j,k,6)
            file_buffer_dble(ix+6) = q_3d(i,j,k,7)
            
            IF(npass .GT. 0) THEN
                DO ipass = 1, npass
                    file_buffer_dble(ix+6+ipass) = qpass_3d(i,j,k,ipass)
                END DO
                ix = ix + npass
            END IF
            
            ix = ix + 7
            
        END DO    
        END DO
        END DO
    
    ELSE

        ! load up output data into the io buffer 
        DO k = 1, nz    
        DO j = 1, ny
        DO i = 1, nx
            
            file_buffer_sngl(ix+0) = SNGL(q_3d(i,j,k,1))   !######WARNING: Need to be more careful with casting DOUBLE to FLOATS. Maybe put in a floor value to avoid floating point exceptions
            file_buffer_sngl(ix+1) = SNGL(q_3d(i,j,k,2))
            file_buffer_sngl(ix+2) = SNGL(q_3d(i,j,k,3))
            file_buffer_sngl(ix+3) = SNGL(q_3d(i,j,k,4))
            file_buffer_sngl(ix+4) = SNGL(q_3d(i,j,k,5))
            file_buffer_sngl(ix+5) = SNGL(q_3d(i,j,k,6))
            file_buffer_sngl(ix+6) = SNGL(q_3d(i,j,k,7))

            IF(npass .GT. 0) THEN
                DO ipass = 1, npass
                    file_buffer_sngl(ix+6+ipass) = SNGL(qpass_3d(i,j,k,ipass))
                END DO
                ix = ix + npass
            END IF
            
            ix = ix + 7
            
        END DO    
        END DO
        END DO
        
    END IF
    
END SUBROUTINE load_parallel_file_buffer


! This subroutine writes the MHD state grid array into a restart file
SUBROUTINE write_restart(t_sim)

    REAL*8, INTENT(IN) :: t_sim
    CHARACTER(LEN=180) :: filename1
    CHARACTER(LEN=6) :: uniti, coordx, coordy, coordz


    IF(myrank .EQ. 0) PRINT*,''
    IF(myrank .EQ. 0) PRINT*,'Writing to restart file...'
    IF(myrank .EQ. 0) PRINT*,''
  
    WRITE(coordx,'(I1.1)') my_coord(1)
    WRITE(coordy,'(I1.1)') my_coord(2)
    WRITE(coordz,'(I1.1)') my_coord(3)
  
 
    filename1 = TRIM(output_filepath)//TRIM('/Restart/restart_mhd_rank=')//TRIM(coordx)//TRIM(coordy) &
                //TRIM(coordz)//TRIM('.dat')
                
    OPEN(UNIT=10,FILE=filename1, FORM = 'UNFORMATTED', ACCESS = 'SEQUENTIAL')
   
    WRITE(10) t_sim
    WRITE(10) q_3d
    WRITE(10) bface_3d
    IF(npass .GT. 0) WRITE(10) qpass_3d
    WRITE(10) dump_count
    WRITE(10) s1   
    WRITE(10) s2   
    WRITE(10) s3   

    CLOSE(UNIT=10)


END SUBROUTINE write_restart


! This subroutine reads initial data from an exisiting restart file
SUBROUTINE read_restart(t_sim)

    REAL*8, INTENT(OUT) :: t_sim
    CHARACTER(LEN=180) :: filename1
    CHARACTER(LEN=6) :: uniti, coordx, coordy, coordz

  
    IF(myrank .EQ. 0) THEN 
        PRINT*,''
        PRINT*,'Reading from restart file...'
        PRINT*,''
    END IF
   
    WRITE(coordx,'(I1.1)') my_coord(1)
    WRITE(coordy,'(I1.1)') my_coord(2)
    WRITE(coordz,'(I1.1)') my_coord(3)
  
 
 
    filename1 = TRIM(output_filepath)//TRIM('/Restart/restart_mhd_rank=')//TRIM(coordx)//TRIM(coordy) &
                //TRIM(coordz)//TRIM('.dat')

    OPEN(UNIT=10,FILE=filename1, FORM = 'UNFORMATTED', STATUS='OLD', ACCESS = 'SEQUENTIAL')

    READ(10) t_sim   
    READ(10) q_3d
    READ(10) bface_3d
    IF(npass .GT. 0) READ(10) qpass_3d
    READ(10) dump_count
    READ(10) s1
    READ(10) s2
    READ(10) s3
    
    CLOSE(UNIT=10)

    PRINT*,''
    PRINT*,'Completed reading from restart file.'
    PRINT*,''
    

END SUBROUTINE read_restart



END MODULE io_mod
