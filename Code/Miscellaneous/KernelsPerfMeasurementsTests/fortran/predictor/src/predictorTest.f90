SUBROUTINE predictorTest(lqhi,lFhi,luh) 
    USE typesDef
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)  :: luh(nVar,nDOF(1),nDOF(2),nDOF(3))              ! spatial degrees of freedom 
    REAL, INTENT(OUT) :: lqhi(nVar,nDOF(1),nDOF(2),nDOF(3))             ! time-averaged space-time degrees of freedom 
    REAL, INTENT(OUT) :: lFhi(nVar,d,nDOF(1),nDOF(2),nDOF(3))           ! time-averaged nonlinear flux tensor in each space-time DOF 

    ! Local variables 
    INTEGER :: i,j,k,l,iVar,iDim, iter 
    REAL    :: rhs0(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))               ! contribution of the initial condition to the known right hand side 
    REAL    :: rhs(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! known right hand side 
    REAL    :: lqh(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! space-time degrees of freedom 
    REAL    :: lFh(nVar,d,nDOF(1),nDOF(2),nDOF(3),nDOF(0))              ! nonlinear flux tensor in each space-time DOF 
    REAL    :: aux(d), w                                                ! auxiliary variables 
    REAL    :: lqhold(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))             ! old space-time degrees of freedom 
    REAL    :: lqx(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! spatial derivative qx of q 
    REAL    :: lqy(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! spatial derivative qy of q 
    REAL    :: lqz(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! spatial derivative qz of q 
    REAL    :: lqt(nVar,nDOF(1),nDOF(2),nDOF(3),nDOF(0))                ! time derivative qt of q 
    REAL    :: res                                                      ! residual 
    REAL, PARAMETER :: tol = 1e-7                                      ! tolerance 

    !
    ! Immediately compute the time-averaged space-time polynomials 
    !
    DO k = 1, nDOF(3)  
     DO j = 1, nDOF(2)  
      DO i = 1, nDOF(1) 
         lqhi(:,i,j,k) = MATMUL( lqh(:,i,j,k,:), wGPN ) 
         DO iDim = 1, nDim 
            lFhi(:,iDim,i,j,k) = MATMUL( lFh(:,iDim,i,j,k,:), wGPN ) 
         ENDDO
      ENDDO
     ENDDO
    ENDDO

    CONTINUE
    !
END SUBROUTINE predictorTest 
    
    

    
    
    
    