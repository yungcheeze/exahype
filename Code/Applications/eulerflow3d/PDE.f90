SUBROUTINE PDEEigenvalues(Lambda,Q,nv)
    USE typesDef, ONLY : nVar, d, EQNgamma 
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE
    ! Argument list 
    REAL, INTENT(IN)  :: Q(nVar), nv(d) 
    REAL, INTENT(OUT) :: Lambda(nVar) 
    ! Local variables 
    REAL :: p, u, c 
    !
    u = ( Q(2)*nv(1) + Q(3)*nv(2) + Q(4)*nv(3) )/Q(1)       ! normal velocity 
    p = (EQNgamma-1)*( Q(5) - 0.5*SUM(Q(2:4)**2)/Q(1) )    ! fluid pressure 
    c = SQRT(EQNgamma*p/Q(1))                              ! sound speed
    !
    Lambda = (/ u-c, u, u, u, u+c /)                        ! The eigenvalues of the Euler equations 
    !
END SUBROUTINE PDEEigenvalues


SUBROUTINE PDEFlux(F,Q)
    USE typesDef, ONLY : nVar, d, EQNgamma 
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE
    ! Argument list 
   REAL, INTENT(IN)  :: Q(nVar) 
    REAL, INTENT(OUT) :: F(nVar,d) 
    ! Local variables 
    REAL :: p, irho  
    !
    ! 3D compressible Euler equations 
    !
    irho = 1.0/Q(1)
    p = (EQNgamma-1)*( Q(5) - 0.5*SUM(Q(2:4)**2)*irho )
    ! 
    F(1,1) = Q(2) 
    F(2,1) = irho*Q(2)*Q(2) + p 
    F(3,1) = irho*Q(2)*Q(3)
    F(4,1) = irho*Q(2)*Q(4)
    F(5,1) = irho*Q(2)*(Q(5)+p)  
    !
    F(1,2) = Q(3) 
    F(2,2) = irho*Q(3)*Q(2)  
    F(3,2) = irho*Q(3)*Q(3) + p 
    F(4,2) = irho*Q(3)*Q(4)
    F(5,2) = irho*Q(3)*(Q(5)+p)  
    ! 
    F(1,3) = Q(4) 
    F(2,3) = irho*Q(4)*Q(2)  
    F(3,3) = irho*Q(4)*Q(3)  
    F(4,3) = irho*Q(4)*Q(4) + p
    F(5,3) = irho*Q(4)*(Q(5)+p)  
     
END SUBROUTINE PDEFlux
    
    
    

 


SUBROUTINE PDENCP(BgradQ,Q,gradQ) 
  USE typesDef, ONLY : nVar, d 
  USE, INTRINSIC :: ISO_C_BINDING 
  IMPLICIT NONE 
  ! Argument list  
  REAL, INTENT(IN)  :: Q(nVar), gradQ(nVar,d) 
  REAL, INTENT(OUT) :: BgradQ(nVar,d) 
  ! Local variables  
  !
  !@todo Please implement
  !
  BgradQ(1, 1) = 0.0
  BgradQ(2, 1) = 0.0
  BgradQ(3, 1) = 0.0
  BgradQ(4, 1) = 0.0
  BgradQ(5, 1) = 0.0
  !
  BgradQ(1, 2) = 0.0
  BgradQ(2, 2) = 0.0
  BgradQ(3, 2) = 0.0
  BgradQ(4, 2) = 0.0
  BgradQ(5, 2) = 0.0
  !
  BgradQ(1, 3) = 0.0
  BgradQ(2, 3) = 0.0
  BgradQ(3, 3) = 0.0
  BgradQ(4, 3) = 0.0
  BgradQ(5, 3) = 0.0
  !
END SUBROUTINE PDENCP 
 


SUBROUTINE PDEMatrixB(Bn,Q,nv) 
  USE typesDef, ONLY : nVar, d 
  USE, INTRINSIC :: ISO_C_BINDING 
  IMPLICIT NONE 
  ! Argument list  
  REAL, INTENT(IN)  :: Q(nVar), nv(d) 
  REAL, INTENT(OUT) :: Bn(nVar,nVar)  
  ! Local variables  
  !
  !@todo Please implement
  !
  Bn(1, 1) = 0.0
  Bn(1, 2) = 0.0
  Bn(1, 3) = 0.0
  Bn(1, 4) = 0.0
  Bn(1, 5) = 0.0
  !
  Bn(2, 1) = 0.0
  Bn(2, 2) = 0.0
  Bn(2, 3) = 0.0
  Bn(2, 4) = 0.0
  Bn(2, 5) = 0.0
  !
  Bn(3, 1) = 0.0
  Bn(3, 2) = 0.0
  Bn(3, 3) = 0.0
  Bn(3, 4) = 0.0
  Bn(3, 5) = 0.0
  !
  Bn(4, 1) = 0.0
  Bn(4, 2) = 0.0
  Bn(4, 3) = 0.0
  Bn(4, 4) = 0.0
  Bn(4, 5) = 0.0
  !
  Bn(5, 1) = 0.0
  Bn(5, 2) = 0.0
  Bn(5, 3) = 0.0
  Bn(5, 4) = 0.0
  Bn(5, 5) = 0.0
  !
  !
END SUBROUTINE PDEMatrixB 
