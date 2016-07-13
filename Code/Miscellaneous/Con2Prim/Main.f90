Program StandaloneC2P
	USE typesDef, ONLY : nVar
	IMPLICIT NONE

	REAL, PARAMETER :: gamma=5./3., t_f=0.4, xsep=0.5
	REAL, PARAMETER :: rho_L= 1, v_L=-0.6, p_L=10
	REAL, PARAMETER :: rho_R=10, v_R=+0.5, p_R=20
	REAL :: Vl0(nVar), Vr0(nVar), Ql(nVar), Qr(nVar), Vl1(nVar), Vr1(nVar)
	INTEGER:: il, ir

	PRINT *, "Playback whats happening when compositing the Shocktube"
	PRINT *, "This program tests whether con2prim is the inverse of prim2con"
	PRINT *, "The 4-vector colums are just V=(rho,vx,vy,vz,p) and Q=(D,sx,sy,sz,e)"

	Vl0(1:nVar) = (/ rho_L, v_L, 0., 0., p_L /)
	Vr0(:) = (/ rho_R, v_R, 0., 0., p_R /)

	CALL PDEPrim2Cons(Ql, Vl0)
	CALL PDEPrim2Cons(Qr, Vr0)

	CALL PDECons2Prim(Vl1, Ql, il)
	CALL PDECons2Prim(Vr1, Qr, il)

	! now of course it should be Vl0==vl1 and Vr0==Vr1.

	PRINT *, "Initial Shocktube primitive variables:"
	PRINT *, "Vl0 =", Vl0
	PRINT *, "Vr0 =", Vr0
	PRINT *, "Intermediate Shocktube conserved variables:"
	PRINT *, "Ql  =", Ql
	PRINT *, "Qr  =", Qr
	PRINT *, "Reconstructed primitive variables:"
	PRINT *, "Vl1 =", Vl1
	PRINT *, "Vr1 =", Vr1

END ! Program StandloneC2P