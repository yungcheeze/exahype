SUBROUTINE MassAccretionRate(Q,masschange)
  USE Parameters, ONLY: gamma, nVar, nDim
  IMPLICIT NONE
  ! Argument list declaration
  REAL :: Q(nVar), masschange
  INTENT(IN)  :: Q
  INTENT(OUT) :: masschange
  ! Local variable declaration
  REAL :: Prim(nVar), Buf(nVar)
  REAL :: rho,vx,vy,vz,p,bx,by,bz,ex,ey,ez,cs,c0
  REAL :: v2,b2,e2,lf,w,ww,uem,gamma1
  REAL :: lapse, gp, gm, dcs, dc0, eel
  REAL :: g_contr(3,3), g_cov(3,3)
  REAl :: shift(3), vf(3), vf_cov(3), Ev(3), Bv(3), ExB(3)
  REAL :: A(3,3), devG(3,3), G(3,3), temp(3,3), Id(3,3), detA, eh, S, evv, T, falpha
  REAL :: alphas, alphal, rhos, rhol, us, vs, ul, vl, es, el
  REAL :: EE(3),BB(3),vv(3),detvEB, vxE(3), vxB(3)
  REAL :: psi, BV_contr(3), Qv_contr(3), QB_contr(3), vxB_contr(3), vb_conv(3), b2_cov, vb_cov
  INTEGER :: i

  ! TODO: Alejandro

  masschange = 0.0
  
END SUBROUTINE MassAccretionRate
