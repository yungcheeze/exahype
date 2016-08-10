PROGRAM Main
    USE typesDef
    IMPLICIT NONE
    
    ! Local variables
    INTEGER :: i,j,k,iVar  
    DOUBLE PRECISION, ALLOCATABLE  :: lduh(:,:,:,:)    ! lduh(nVar,nDOF(1),nDOF(2),nDOF(3))
    DOUBLE PRECISION, ALLOCATABLE  :: lFbnd(:,:,:,:)   ! time-averaged nonlinear flux tensor in each space-time DOF, lFbnd(nVar,6,nDOF(2),nDOF(3))
    
    ! Input data
    CHARACTER(len=64) :: arg
    CALL GET_COMMAND_ARGUMENT(1, arg)
    OPEN (unit=99, file=arg, form='formatted', status='old', action='read')

    ! Gauss-Legendre weights
    CALL ADERDGInit
    
    ! initialize input data
    ALLOCATE(lFbnd(nVar,6, nDOF(2),nDOF(3)))
    ALLOCATE(lduh(nVar,nDOF(1),nDOF(2),nDOF(3)))
    
    DO k=1,nDOF(3)
      DO j=1,nDOF(2)
        DO i=1,nDOF(1)
          READ(99,*) lduh(:,i,j,k)
        ENDDO
      ENDDO
    ENDDO
   
    DO k=1,nDOF(3)
      DO j=1,nDOF(2)
        DO iVar=1,nVar
          READ(99,*) lFbnd(iVar,:,j,k)
        ENDDO
      ENDDO
    ENDDO
   
    !run and output
    CALL ADERSurfaceIntegral(lduh, lFbnd, dx)
    PRINT *, lduh
 
    DEALLOCATE(lFbnd)
    DEALLOCATE(lduh)

END PROGRAM Main
