PROGRAM Main
    USE typesDef
    IMPLICIT NONE
    ! Local variables
    INTEGER :: i,j,k,l,iElem,iFace,iVar,di
    REAL, POINTER :: lQbndL(:,:,:),lFbndL(:,:,:),lQbndR(:,:,:),lFbndR(:,:,:)
    REAL :: cnt,cntk,cntj,cnti
    ! Input data
    CHARACTER(len=64) :: arg
    CALL GET_COMMAND_ARGUMENT(1, arg)
    OPEN (unit=99, file=arg, form='formatted', status='old', action='read')
    !
    ! Init variables for predictor

    ! Gauss-Legendre weights
    CALL ADERDGInit
    
    !lFhi(nVar,d,nDOF(1),nDOF(2),nDOF(3))
    
    ! initialize input data

!    cnt = 1.d0   
!
!      DO k= 1, nDOF(3) ! 1 -> loop redundant
!        DO j = 1,nDOF(2)
!          DO i=1, nDOF(1)
!            DO iVar = 1,nVar
!              Fhi(iVar,1,i,j,k,1) = cnt
!              cnt = cnt+1.d0
!            ENDDO
!          ENDDO
!        ENDDO
!      ENDDO
!DO k= 1, nDOF(3) ! 1 -> loop redundant
!        DO i = 1,nDOF(1)
!          DO j=1,nDOF(2)
!            DO iVar = 1,nVar
!              Fhi(iVar,2,i,j,k,1) = cnt
!              cnt = cnt+1.d0
!            ENDDO
!          ENDDO
!        ENDDO
!      ENDDO
    
    DO di =1,2
      DO iVar = 1,nVar
        DO i = 1,nDOF(1)
          READ(99,*) Fhi(iVar,di,i,:,1,1)
        ENDDO
      ENDDO    
    ENDDO
    
      
   
    !PRINT *, Fhi(1,1,:,:,1,1)
    !PRINT *, Fhi(3,2,3,4,1,1)
   
   
    CALL ADERVolumeIntegral(duh(:,:,:,:,1),qhi(:,:,:,:,1),Fhi(:,:,:,:,:,1))  
    PRINT *, duh(:,:,:,1,1)

END PROGRAM Main
