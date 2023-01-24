MODULE readfile_mod

USE constants_mod
IMPLICIT NONE

CONTAINS



SUBROUTINE readfile_ranks_multivar(t_dump, npass, varlow, varhi, double_prec, field_array)


    INTEGER, INTENT(IN) :: t_dump, npass, varlow, varhi
	LOGICAL, INTENT(IN) :: double_prec
    REAL(4) :: field_array(1:nx*nranks_x,1:ny*nranks_y,1:nz*nranks_z,1:1+varhi-varlow) ! I'm making the field array single precision by default
    INTEGER :: i, j, k, ix, iy, iz, iv, rankx, ranky, rankz, ilow, ihi, jlow, jhi, klow, khi
    INTEGER :: nvars, item_bytesize, byte_offset
    CHARACTER(LEN=300) :: filename
    CHARACTER(LEN=6) :: uniti, rx, ry, rz
    CHARACTER(LEN=20) :: riemann_solver_prefix
    REAL(8) :: temp
    
    IF(double_prec) THEN
        item_bytesize = 8
    ELSE
        item_bytesize = 4
    END IF   
    
    nvars = 7 + npass
    
    PRINT*,'Reading file for dump # ',t_dump
    PRINT*,''
   
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
      
    
    DO rankz = 0, nranks_z-1 
    DO ranky = 0, nranks_y-1
    DO rankx = 0, nranks_x-1
    
        WRITE(rx,'(I1.1)') rankx
        WRITE(ry,'(I1.1)') ranky
        WRITE(rz,'(I1.1)') rankz

        filename = TRIM(output_filepath)//TRIM('/Snapshots/')//TRIM('xyzcube_rank=')//TRIM(rx)//TRIM(ry)// & 
                   TRIM(rz)//TRIM('_dump=')//TRIM(uniti)//TRIM('.dat')
                
        !PRINT*,''        
        !PRINT*,'Reading from file : '
        !WRITE(*,*) filename        
        
        OPEN(UNIT=17,FILE=filename, FORM = 'UNFORMATTED', STATUS = 'OLD', ACCESS = 'STREAM')
    
        
        ilow = 1 + rankx * nx 
        ihi  = ilow + nx - 1
        jlow = 1 + ranky * ny
        jhi  = jlow + ny - 1
        klow = 1 + rankz * nz
        khi  = klow + nz - 1
        
        DO k = klow, khi 
        iz = k - rankz * nz
        DO j = jlow, jhi
        iy = j - ranky * ny
        DO i = ilow, ihi
        ix = i - rankx * nx

            byte_offset =  1 + (varlow-1)*item_bytesize +  (ix-1) * (nvars*item_bytesize) + (iy-1) * nx * (nvars*item_bytesize) + (iz-1) * nx * ny * (nvars*item_bytesize) 

            DO iv = varlow, varhi
           
                IF(.NOT. double_prec) THEN   
                   READ(17, POS = byte_offset) field_array(i,j,k,iv-varlow+1)
                ELSE
                   READ(17, POS = byte_offset) temp
                   field_array(i,j,k,iv-varlow+1) = SNGL(temp)
                END IF
            
                byte_offset = byte_offset + item_bytesize
           
            END DO
           
               
        END DO    
        END DO
        END DO

        CLOSE(UNIT=17)
    
    END DO    
    END DO
    END DO    
    
    
    
    
END SUBROUTINE readfile_ranks_multivar


SUBROUTINE readfile_ranks_singlevar(t_dump, npass, varpos, double_prec, field_array)


    INTEGER, INTENT(IN) :: t_dump, npass, varpos
	LOGICAL, INTENT(IN) :: double_prec
    REAL(4), INTENT(INOUT) :: field_array(1:nx*nranks_x,1:ny*nranks_y,1:nz*nranks_z) ! I'm making the field array single precision by default
    INTEGER :: i, j, k, ix, iy, iz,  rankx, ranky, rankz, ilow, ihi, jlow, jhi, klow, khi
    INTEGER :: nvars, item_bytesize, byte_offset
    CHARACTER(LEN=300) :: filename
    CHARACTER(LEN=6) :: uniti, rx, ry, rz
    CHARACTER(LEN=20) :: riemann_solver_prefix
    REAL(8) :: temp

    
    IF(double_prec) THEN
        item_bytesize = 8
    ELSE
        item_bytesize = 4
    END IF   
    
    nvars = 7 + npass
    

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
      
    
    
    DO rankz = 0, nranks_z-1 
    DO ranky = 0, nranks_y-1
    DO rankx = 0, nranks_x-1
    
        WRITE(rx,'(I1.1)') rankx
        WRITE(ry,'(I1.1)') ranky
        WRITE(rz,'(I1.1)') rankz

        filename = TRIM(output_filepath)//TRIM('/Snapshots/')//TRIM('xyzcube_rank=')//TRIM(rx)//TRIM(ry)// & 
                   TRIM(rz)//TRIM('_dump=')//TRIM(uniti)//TRIM('.dat')
                
        !PRINT*,''        
        !PRINT*,'Reading from file : '
        !WRITE(*,*) filename        
                
        OPEN(UNIT=17,FILE=filename, FORM = 'UNFORMATTED', STATUS = 'OLD', ACCESS = 'STREAM')
    
        
        ilow = 1 + rankx * nx 
        ihi  = ilow + nx - 1
        jlow = 1 + ranky * ny
        jhi  = jlow + ny - 1
        klow = 1 + rankz * nz
        khi  = klow + nz - 1
    
    
        DO k = klow, khi 
        iz = k - rankz * nz
        DO j = jlow, jhi
        iy = j - ranky * ny
        DO i = ilow, ihi
        ix = i - rankx * nx

            byte_offset =  1 + (varpos-1)*item_bytesize +  (ix-1) * (nvars*item_bytesize) + (iy-1) * nx * (nvars*item_bytesize) + (iz-1) * nx * ny * (nvars*item_bytesize) 

            IF(.NOT. double_prec) THEN   
               READ(17, POS = byte_offset) field_array(i,j,k)
            ELSE
               READ(17, POS = byte_offset) temp
               field_array(i,j,k) = SNGL(temp)
            END IF
           
           
        END DO    
        END DO
        END DO

        CLOSE(UNIT=17)
    
    END DO    
    END DO
    END DO    
    
    
    
    
END SUBROUTINE readfile_ranks_singlevar







END MODULE readfile_mod 