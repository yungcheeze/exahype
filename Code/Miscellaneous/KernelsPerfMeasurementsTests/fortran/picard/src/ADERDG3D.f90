PROGRAM ADERDG3D 
    USE typesDef
    IMPLICIT NONE
    ! Local variables
    INTEGER :: i,j,k,iElem,iFace
    REAL, POINTER :: lQbndL(:,:,:),lFbndL(:,:,:),lQbndR(:,:,:),lFbndR(:,:,:)

    ! We first need to compute the relevant matrices, set initial
    ! conditions and prepare some necessary stuff...  
    CALL ADERDGInit 
    
    CALL CPU_TIME(tCPU1) 
    ! Main loop in time 
    DO timestep = 1, 1000
      CALL ADERSpaceTimePredictor(qhi(:,:,:,:,1),Fhi(:,:,:,:,:,1),uh(:,:,:,:,1))   
    ENDDO    
    CALL CPU_TIME(tCPU2)
    
    TEU = timestep*nElem 
    PRINT *, N,' ', tCPU2-tCPU1 
    
    
END PROGRAM ADERDG3D
    