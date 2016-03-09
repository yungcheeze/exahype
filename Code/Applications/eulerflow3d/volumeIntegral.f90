SUBROUTINE ADERVolumeIntegral(lduh,lFhi,dx)
    USE typesDef

    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)              :: lFhi(nVar,nDOF(1),nDOF(2),nDOF(3),d)    ! nonlinear flux tensor in each space-time DOF 
    REAL, INTENT(OUT)             :: lduh(nVar,nDOF(1),nDOF(2),nDOF(3))      ! spatial degrees of freedom 
    DOUBLE PRECISION, INTENT(IN)  :: dx(d)                                          ! 
    ! Local variables 
    INTEGER           :: i,j,k,l
    REAL              :: aux(d) 
    ! 
    ! Initialize the update DOF 
    lduh = 0. 
    ! x - direction 

    !PRINT *, ' --------lFhi-volumeIntegral-------------------------------- ' 
    !PRINT *, lFhi
    !PRINT *, ' --------lFhi-volumeIntegral-------------------------------- ' 
    !CALL EXIT    

    DO k = 1, nDOF(3)
        DO j = 1, nDOF(2) 
            aux = (/ 1.d0, wGPN(j), wGPN(k) /) 
            lduh(:,:,j,k) = lduh(:,:,j,k) + MATMUL( lFhi(:,:,j,k,1), TRANSPOSE(Kxi) )*PRODUCT(aux(1:nDim))/dx(1) 
        ENDDO
    ENDDO
    
    ! y - direction (not needed in 1D) 
    IF(nDim>=2) THEN
        DO k = 1, nDOF(3)
            DO i = 1, nDOF(1) 
                aux = (/ 1.d0, wGPN(i), wGPN(k) /) 
                lduh(:,i,:,k) = lduh(:,i,:,k) + MATMUL( lFhi(:,i,:,k,2), TRANSPOSE(Kxi) )*PRODUCT(aux(1:nDim))/dx(2)
            ENDDO
        ENDDO
    ENDIF 
    
    ! z - direction (node needed in 1D and 2D) 
    IF(nDim>=3) THEN
        DO j = 1, nDOF(2)
            DO i = 1, nDOF(1)
                aux = (/ 1.d0, wGPN(i), wGPN(j) /)  
                lduh(:,i,j,:) = lduh(:,i,j,:) + MATMUL( lFhi(:,i,j,:,3), TRANSPOSE(Kxi) )*PRODUCT(aux(1:nDim))/dx(3)  
            ENDDO
        ENDDO
    ENDIF 
    !
    CONTINUE
    !
    !PRINT *, ' -----lduh------------------------------------ ' 
    !PRINT *, lduh
    !PRINT *, ' -----lduh------------------------------------ ' 

    !OPEN(UNIT=12, FILE="aoutput.txt", ACTION="write", STATUS="replace")
    !WRITE(12, '(ES24.16,1x)') , lduh
    !CALL EXIT    
    
END SUBROUTINE ADERVolumeIntegral 
    
    
