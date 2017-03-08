SUBROUTINE MassAccretionRate(Q,masschange)
  USE Parameters, ONLY: gamma, nVar, nDim
  IMPLICIT NONE
  ! Argument list declaration
  REAL :: Q(nVar), masschange
  INTENT(IN)  :: Q
  INTENT(OUT) :: masschange

  ! Local variable declaration
  REAL :: v2,b2,e2,lf,w,ww,uem,gamma1
  REAL :: lapse, g_contr(3,3), g_cov(3,3), gp
  REAl :: shift(3), vf(3), vf_cov(3), ur_contr(3),V(nVar)
  INTEGER :: i, iErr

  ! TODO: Alejandro

  CALL PDECons2Prim(V,Q,iErr)

  !gammaij = V(14:19)
  g_cov(1,1) = V(14)
  g_cov(1,2) = V(15)
  g_cov(1,3) = V(16)
  g_cov(2,2) = V(17)
  g_cov(2,3) = V(18)
  g_cov(3,3) = V(19)
  g_cov(2,1) = V(15)
  g_cov(3,1) = V(16)
  g_cov(3,2) = V(18)

  ! Compute the inverse of the 3-metric
  CALL MatrixInverse3x3(g_cov,g_contr,gp)

  vf_cov = V(2:4)
  ! v_i --> v^i From covariant to contravariant
  vf     = MATMUL(g_contr,vf_cov)
  gp = SQRT(gp)
  lapse = V(10)
  shift = V(11:13)

  ! u^i = v^i - beta^i/alpha
  !thinking in spherical simmetry
  
  ur_contr = vf(1)  - shift(1)/lapse
  masschange = lapse *gp* Q(1) * ur_contr(1)
  
  
END SUBROUTINE MassAccretionRate
