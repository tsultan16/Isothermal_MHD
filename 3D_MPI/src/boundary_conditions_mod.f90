!TODO: Right now, I only included all boundaries either open or periodic. Need to allow 
!      choice of open/periodic/reflecting for each boundary separately. 
!

MODULE boundary_conditions_mod


USE constants_mod
USE grid_data_mod

IMPLICIT NONE


CONTAINS


! Apply global boundary conditions to state vector 
! (Note: boundary conditions are only appled to the edges of the local MPI domain that
! coincide with the global computational domain.)
SUBROUTINE boundary_conditions()

     INTEGER :: i, j, k, nb_pass
     
     nb_pass = 3
     
     ! Open Boundaries
     IF(boundary_type .EQ. 1) THEN
    
        !************        
        ! x boundary
        !************        
        
        IF(my_coord(1) .EQ. 0) THEN
        
            DO i = 0, nb-1
                q_3d(1-nb+i,:,:,:) = q_3d(1,:,:,:)          
                bface_3d(1-nb+i,:,:,:) = bface_3d(1,:,:,:)         
            END DO

            IF(npass .GT. 0) THEN
                DO i = 0, nb_pass-1
                    qpass_3d(1-nb_pass+i,:,:,:) = qpass_3d(1,:,:,:)
                END DO
            END IF
            
        END IF
        
        IF(my_coord(1) .EQ. nranks_x-1) THEN
        
            DO i = 0, nb-1
                q_3d(nx+1+i,:,:,:) = q_3d(nx,:,:,:)      
                bface_3d(nx+1+i,:,:,:) = bface_3d(nx,:,:,:)           
            END DO

            IF(npass .GT. 0) THEN
                DO i = 0, nb_pass-1
                    qpass_3d(nx+1+i,:,:,:) = qpass_3d(nx,:,:,:)
                END DO
            END IF
                
        END IF
        
        
        !************
        ! y boundary
        !************

        IF(my_coord(2) .EQ. 0) THEN

            DO j = 0, nb-1
                q_3d(:,1-nb+j,:,:) = q_3d(:,1,:,:)                   
                bface_3d(:,1-nb+j,:,:) = bface_3d(:,1,:,:)     
            END DO   
          
            IF(npass .GT. 0) THEN
                DO j = 0, nb_pass-1
                    qpass_3d(:,1-nb_pass+j,:,:) = qpass_3d(:,1,:,:) 
                END DO
            END IF
            
        END IF

        IF(my_coord(2) .EQ. nranks_y-1) THEN

            DO j = 0, nb-1
                q_3d(:,ny+1+j,:,:) = q_3d(:,ny,:,:)                     
                bface_3d(:,ny+1+j,:,:) = bface_3d(:,ny,:,:)   
            END DO      
        
            IF(npass .GT. 0) THEN
                DO j = 0, nb_pass-1
                    qpass_3d(:,ny+1+j,:,:) = qpass_3d(:,ny,:,:) 
                END DO
            END IF
        
        END IF


        !************ 
        ! z boundary
        !************

        IF(my_coord(3) .EQ. 0) THEN

            DO k = 0, nb-1
                q_3d(:,:,1-nb+k,:) = q_3d(:,:,1,:)         
                bface_3d(:,:,1-nb+k,:) = bface_3d(:,:,1,:)     
            END DO

            IF(npass .GT. 0) THEN
                DO k = 0, nb_pass-1
                    qpass_3d(:,:,1-nb_pass+k,:) = qpass_3d(:,:,1,:)    
                END DO
            END IF            
            
        END IF
        
        IF(my_coord(3) .EQ. nranks_z-1) THEN

            DO k = 0, nb-1
                q_3d(:,:,nz+1+k,:) = q_3d(:,:,nz,:)      
                bface_3d(:,:,nz+1+k,:) = bface_3d(:,:,nz,:)      
            END DO
        
            IF(npass .GT. 0) THEN
                DO k = 0, nb_pass-1
                    qpass_3d(:,:,nz+1+k,:) = qpass_3d(:,:,nz,:)      
                END DO
            END IF 
                
        END IF
         
    END IF    


    ! Periodic Boundaries
    IF(boundary_type .EQ. 2) THEN
    
        !************        
        ! x boundary
        !************        
        
        IF(nranks_x .EQ. 1) THEN
        
            DO i = 0, nb-1        
                q_3d(-i,:,:,:) = q_3d(nx-i,:,:,:)          
                bface_3d(-i,:,:,:) = bface_3d(nx-i,:,:,:)         
            END DO

            IF(npass .GT. 0) THEN
                DO i = 0, nb_pass-1
                    qpass_3d(-i,:,:,:) = qpass_3d(nx-i,:,:,:)     
                END DO
            END IF    

            DO i = 0, nb-1       
                q_3d(nx+1+i,:,:,:) = q_3d(1+i,:,:,:)      
                bface_3d(nx+1+i,:,:,:) = bface_3d(1+i,:,:,:)           
            END DO

            IF(npass .GT. 0) THEN
                DO i = 0, nb_pass-1
                    qpass_3d(nx+1+i,:,:,:) = qpass_3d(1+i,:,:,:)      
                END DO
            END IF    
                
        END IF
        
        
        !************
        ! y boundary
        !************

        IF(nranks_y .EQ. 1) THEN

            DO j = 0, nb-1
                q_3d(:,-j,:,:) = q_3d(:,ny-j,:,:)                   
                bface_3d(:,-j,:,:) = bface_3d(:,ny-j,:,:)     
            END DO   

            IF(npass .GT. 0) THEN
                DO j = 0, nb_pass-1
                    qpass_3d(:,-j,:,:) = qpass_3d(:,ny-j,:,:)        
                END DO
            END IF 
            
            DO j = 0, nb-1
                q_3d(:,ny+1+j,:,:) = q_3d(:,1+j,:,:)                     
                bface_3d(:,ny+1+j,:,:) = bface_3d(:,1+j,:,:)    
            END DO      
        
            IF(npass .GT. 0) THEN
                DO j = 0, nb_pass-1
                    qpass_3d(:,ny+1+j,:,:) = qpass_3d(:,1+j,:,:)        
                END DO
            END IF 
            
        END IF


        !************ 
        ! z boundary
        !************

        IF(nranks_z .EQ. 1) THEN

            DO k = 0, nb-1
                q_3d(:,:,-k,:) = q_3d(:,:,nz-k,:)         
                bface_3d(:,:,-k,:) = bface_3d(:,:,nz-k,:)     
            END DO

            IF(npass .GT. 0) THEN
                DO k = 0, nb_pass-1
                    qpass_3d(:,:,-k,:) = qpass_3d(:,:,nz-k,:)       
                END DO
            END IF 
            
            DO k = 0, nb-1
                q_3d(:,:,nz+1+k,:) = q_3d(:,:,1+k,:)      
                bface_3d(:,:,nz+1+k,:) = bface_3d(:,:,1+k,:)      
            END DO
        
            IF(npass .GT. 0) THEN
                DO k = 0, nb_pass-1
                    qpass_3d(:,:,nz+1+k,:) = qpass_3d(:,:,1+k,:)   
                END DO
            END IF 
        
        END IF            

    END IF
    
        
END SUBROUTINE boundary_conditions



END MODULE boundary_conditions_mod
