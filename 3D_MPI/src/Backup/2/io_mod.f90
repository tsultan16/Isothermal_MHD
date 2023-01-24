MODULE io_mod

USE constants_mod
USE grid_data_mod
USE mpi_domain_mod

IMPLICIT NONE

CONTAINS


! This subroutine writes the MHD state variables grid array from a single time-step 
! into an unformatted (i.e. binary) file, i.e. a stream of bytes.
! Each MPI rank writes to it's on file.
SUBROUTINE writetofile_unformatted(tstep)


    INTEGER, INTENT(IN) :: tstep
	
    INTEGER :: i, j, k
    REAL*8 :: dx, vx, vy, vz, vpar, vperp, bpar, bperp, isqr2
    REAL*8 :: x, y, z
    CHARACTER(LEN=300) :: filename1, filename2, filename3, filename4, filename5, filename6, filename7
    CHARACTER(LEN=6) :: uniti, coordx, coordy, coordz
    
    
    IF(tstep<10) THEN
        WRITE(uniti,'(I1.1)') tstep
    ELSE IF(tstep>=10 .and. tstep<100) THEN
        WRITE(uniti,'(I2.2)') tstep
    ELSE IF(tstep>=100 .and. tstep<1000) THEN
        WRITE (uniti,'(I3.3)') tstep
    ELSE IF(tstep>=1000 .and. tstep<10000) THEN
        WRITE (uniti,'(I4.3)') tstep
    ELSE IF(tstep>=10000 .and. tstep<100000) THEN
        WRITE (uniti,'(I5.3)') tstep  
    END IF
      
    
    WRITE(coordx,'(I1.1)') my_coord(1)
    WRITE(coordy,'(I1.1)') my_coord(2)
    WRITE(coordz,'(I1.1)') my_coord(3)
    
       
    filename1 = TRIM(output_filepath)//TRIM('xcut_rank=')//TRIM(coordx)//TRIM(coordy) &
                //TRIM(coordz)//TRIM('_t=')//TRIM(uniti)//TRIM('.dat')
                
    filename2 = TRIM(output_filepath)//TRIM('ycut_rank=')//TRIM(coordx)//TRIM(coordy) &
                //TRIM(coordz)//TRIM('_t=')//TRIM(uniti)//TRIM('.dat')
                
    filename3 = TRIM(output_filepath)//TRIM('zcut_rank=')//TRIM(coordx)//TRIM(coordy) & 
                //TRIM(coordz)//TRIM('_t=')//TRIM(uniti)//TRIM('.dat')
                
    filename4 = TRIM(output_filepath)//TRIM('xydiagcut_rank=')//TRIM(coordx)//TRIM(coordy) &
                //TRIM(coordz)//TRIM('_t=')//TRIM(uniti)//TRIM('.dat')
                
    filename5 = TRIM(output_filepath)//TRIM('yzdiagcut_rank=')//TRIM(coordx)//TRIM(coordy) & 
                //TRIM(coordz)//TRIM('_t=')//TRIM(uniti)//TRIM('.dat')

    filename6 = TRIM(output_filepath)//TRIM('xyplane_rank=')//TRIM(coordx)//TRIM(coordy) & 
                //TRIM(coordz)//TRIM('_t=')//TRIM(uniti)//TRIM('.dat')
 
    filename7 = TRIM(output_filepath)//TRIM('xyzcube_rank=')//TRIM(coordx)//TRIM(coordy) & 
                //TRIM(coordz)//TRIM('_t=')//TRIM(uniti)//TRIM('.dat')
                
    dx = 1.d0 / MAX(nx*nranks_x,ny*nranks_y,nz*nranks_z)
    isqr2 = 1.d0 / SQRT(2.d0)
 

    IF(output_xcut) THEN
 
    OPEN(UNIT=11,FILE=filename1, FORM = 'UNFORMATTED', ACCESS = 'STREAM')
    
    
    j = ny-2
    k = nz/2
    
    DO i =1,nx
        vx = q_3d(i,j,k,2)/q_3d(i,j,k,1)
        vy = q_3d(i,j,k,3)/q_3d(i,j,k,1)
        vz = q_3d(i,j,k,4)/q_3d(i,j,k,1)
        WRITE(11) (my_coord(1)*nx+(i-0.5))*dx,q_3d(i,j,k,1),vx,vy,vz,q_3d(i,j,k,5),q_3d(i,j,k,6),q_3d(i,j,k,7)
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
        WRITE(12) my_coord(2)+(j-0.5)*dx,q_3d(i,j,k,1),vx,vy,vz,q_3d(i,j,k,5),q_3d(i,j,k,6),q_3d(i,j,k,7)
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
        WRITE(13) my_coord(3)+(k-0.5)*dx,q_3d(i,j,k,1),vx,vy,vz,q_3d(i,j,k,5),q_3d(i,j,k,6),q_3d(i,j,k,7)
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
        
            WRITE(14) (my_coord(1)*nx+(i-0.5))*dx,q_3d(i,j,k,1),vperp,vpar,vz,bperp,bpar,q_3d(i,j,k,7)
            
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
        
            WRITE(15) j*dx,q_3d(i,j,k,1),vperp,vpar,vx,bperp,bpar,q_3d(i,j,k,5)
            
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
        WRITE(16) (my_coord(1)*nx+(i-0.5))*dx,(my_coord(2)*ny+(j-0.5))*dx,q_3d(i,j,k,1),vx,vy,vz,q_3d(i,j,k,5),q_3d(i,j,k,6),q_3d(i,j,k,7)
        
    END DO
    END DO
    
    CLOSE(UNIT=16)
    
    END IF
    

    IF(output_xyz_cube) THEN

    OPEN(UNIT=17,FILE=filename7, FORM = 'UNFORMATTED', ACCESS = 'STREAM')
    
    DO k =1,nz 
    DO j =1,ny
    DO i =1,nx
        
        vx = q_3d(i,j,k,2)/q_3d(i,j,k,1)
        vy = q_3d(i,j,k,3)/q_3d(i,j,k,1)
        vz = q_3d(i,j,k,4)/q_3d(i,j,k,1)
        WRITE(17) (my_coord(1)*nx+(i-0.5))*dx,(my_coord(2)*ny+(j-0.5))*dx,q_3d(i,j,k,1),vx,vy,vz,q_3d(i,j,k,5),q_3d(i,j,k,6),q_3d(i,j,k,7)
        
    END DO    
    END DO
    END DO
    
    CLOSE(UNIT=17)
    
    END IF
    
    
END SUBROUTINE writetofile_unformatted


! Parallel version of file output routine. Each mpi rank writes
! to a different portion of a common file visible to all ranks.
! Use this only when a high-performance parallel filesystem is 
! available (e.g. Lustre, GPFS, etc.), otherwise this will be slower
! than each rank writing to it's own file (i.e. don't use this routine 
! if running on a single node).
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
SUBROUTINE writetofile_unformatted_parallel(tstep, parallel_io_req, file_handle, initial, final)


    INTEGER, INTENT(IN) :: tstep
    INTEGER, INTENT(INOUT) :: parallel_io_req, file_handle
	LOGICAL, INTENT(IN) :: initial, final
    INTEGER :: ierr, status(MPI_STATUS_SIZE)
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset
    CHARACTER(LEN=300) :: filename
    CHARACTER(LEN=6) :: uniti

        
    IF(tstep<10) THEN
        WRITE(uniti,'(I1.1)') tstep
    ELSE IF(tstep>=10 .and. tstep<100) THEN
        WRITE(uniti,'(I2.2)') tstep
    ELSE IF(tstep>=100 .and. tstep<1000) THEN
        WRITE (uniti,'(I3.3)') tstep
    ELSE IF(tstep>=1000 .and. tstep<10000) THEN
        WRITE (uniti,'(I4.3)') tstep
    ELSE IF(tstep>=10000 .and. tstep<100000) THEN
        WRITE (uniti,'(I5.3)') tstep  
    END IF
      
    ! wait for the previous file write to complete, then close the old file  
    IF(.NOT. initial) THEN
        CALL MPI_WAIT(parallel_io_req, status, ierr)
        CALL MPI_FILE_CLOSE(file_handle, ierr)
    END IF    
      
    ! now start the next non-blocking file write
    IF(.NOT. final) THEN  
        
        filename = TRIM(output_filepath)//TRIM('mhd_state_common_t=')//TRIM(uniti)//TRIM('.dat')    
        
        ! set the file offset (i.e. number of bytes to skip from the start of the file) 
        offset = ( my_coord(1) + nranks_x * my_coord(2) + nranks_x * nranks_y * my_coord(3) ) * parallel_filesize * 8  ! 8 bytes per data item

        ! open new parallel file
        CALL MPI_FILE_OPEN(comm3d, filename, MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, file_handle, ierr)

        ! write to parallel file (Note: this is a non-blocking MPI operation)
        CALL MPI_FILE_IWRITE_AT(file_handle, offset, file_buffer, parallel_filesize, MPI_DOUBLE_PRECISION, parallel_io_req, ierr)

        IF(ierr .NE. MPI_SUCCESS) THEN
            PRINT*,'Parallel file write failed!'
            !STOP
        END IF

    END IF

END SUBROUTINE writetofile_unformatted_parallel


! this subroutine copies mhd state variables into our parallel file io buffer
SUBROUTINE load_parallel_file_buffer()
	
    INTEGER :: i, j, k, ix
    REAL*8 :: x, y, z, dx
 
    dx = 1.d0 / MAX(nx*nranks_x,ny*nranks_y,nz*nranks_z)
  
    ix = 1

    ! load up output data into the io buffer 
    DO k = 1, nz    
    DO j = 1, ny
    DO i = 1, nx
        
        file_buffer(ix)   = (my_coord(1)*nx+(i-0.5))*dx
        file_buffer(ix+1) = (my_coord(2)*ny+(j-0.5))*dx
        file_buffer(ix+2) = (my_coord(3)*nz+(k-0.5))*dx
        file_buffer(ix+3) = q_3d(i,j,k,1)
        file_buffer(ix+4) = q_3d(i,j,k,2)
        file_buffer(ix+5) = q_3d(i,j,k,3)
        file_buffer(ix+6) = q_3d(i,j,k,4)
        file_buffer(ix+7) = q_3d(i,j,k,5)
        file_buffer(ix+8) = q_3d(i,j,k,6)
        file_buffer(ix+9) = q_3d(i,j,k,7)

        ix = ix  + 10
        
    END DO    
    END DO
    END DO
    
    

END SUBROUTINE load_parallel_file_buffer


! This subroutine writes the MHD state grid array into a restart file
SUBROUTINE write_restart(tstep)

    INTEGER, INTENT(IN) :: tstep
    CHARACTER(LEN=180) :: filename1
    CHARACTER(LEN=6) :: uniti, coordx, coordy, coordz


    PRINT*,''
    IF(myrank .EQ. 0) PRINT*,'Writing to restart file...'
    PRINT*,''
  
    WRITE(coordx,'(I1.1)') my_coord(1)
    WRITE(coordy,'(I1.1)') my_coord(2)
    WRITE(coordz,'(I1.1)') my_coord(3)
  
 
    filename1 = trim('Output/Restart/restart_mhd_rank=')//TRIM(coordx)//TRIM(coordy) &
                //TRIM(coordz)//TRIM('.dat')
                
    OPEN(UNIT=10,FILE=filename1, FORM = 'UNFORMATTED', ACCESS = 'SEQUENTIAL')
   
    WRITE(10) tstep
    WRITE(10) q_3d
    WRITE(10) bface_3d   

    CLOSE(UNIT=10)


END SUBROUTINE write_restart


! This subroutine reads initial data from an exisiting restart file
SUBROUTINE read_restart(tstart)

    INTEGER, INTENT(INOUT) :: tstart
    CHARACTER(LEN=180) :: filename1
    CHARACTER(LEN=6) :: uniti, coordx, coordy, coordz

  
    PRINT*,''
    PRINT*,'Reading from restart file...'
    PRINT*,''
    
    WRITE(coordx,'(I1.1)') my_coord(1)
    WRITE(coordy,'(I1.1)') my_coord(2)
    WRITE(coordz,'(I1.1)') my_coord(3)
  
 
 
    filename1 = trim('Output/Restart/restart_mhd_rank=')//TRIM(coordx)//TRIM(coordy) &
                //TRIM(coordz)//TRIM('.dat')

    OPEN(UNIT=10,FILE=filename1, FORM = 'UNFORMATTED', STATUS='OLD', ACCESS = 'SEQUENTIAL')

    READ(10) tstart    
    READ(10) q_3d
    READ(10) bface_3d
    
    CLOSE(UNIT=10)

    PRINT*,''
    PRINT*,'Completed reading from restart file.'
    PRINT*,''
    

END SUBROUTINE read_restart



END MODULE io_mod
