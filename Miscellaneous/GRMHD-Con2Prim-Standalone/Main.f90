RECURSIVE SUBROUTINE DIFFID (xc, t, gradQ)
  IMPLICIT NONE 
  !
  INTEGER, PARAMETER :: nVar = 19
  REAL, intent(IN)  :: xc(3), t
  REAL, intent(OUT) :: gradQ(3,nVar)
  REAL :: xp(3), xm(3), xp1(3), xm1(3), xp2(3), xm2(3), Qp(nVar), Qm(nVar), Qp1(nVar), Qm1(nVar), Qp2(nVar), Qm2(nVar), V(nVar)
  REAL, PARAMETER   :: epsilon = 1e-7, eps4 = 1e-4
  INTEGER :: i
  !
  ! Metric derivative computed with a fourth order central finite difference 
  !
  DO i = 1, 3
    xp1 = xc
    xp1(i) = xp1(i)+eps4 
    xm1 = xc
    xm1(i) = xm1(i)-eps4 
    xp2 = xc
    xp2(i) = xp2(i)+2*eps4 
    xm2 = xc
    xm2(i) = xm2(i)-2*eps4 
    CALL AlfenWavePrim ( xp1, t, V)
    CALL PDEPrim2Cons(   Qp1, V)
    CALL AlfenWavePrim ( xm1, t, V)
    CALL PDEPrim2Cons(   Qm1, V)
    CALL AlfenWavePrim ( xp2, t, V)
    CALL PDEPrim2Cons(   Qp2, V)
    CALL AlfenWavePrim ( xm2, t, V)
    CALL PDEPrim2Cons(   Qm2, V)

    gradQ(i,:) = ( 8.0*Qp1      -8.0*Qm1       +Qm2       -Qp2       )/(12.0*eps4) 
  ENDDO
END SUBROUTINE DIFFID

Program GRMHDStandaloneC2P
	USE Parameters, ONLY : nVar, nDim, gamma
	IMPLICIT NONE

	REAL :: V0(nVar), V1(nVar), V2(nVar), Q0(nVar), Q1(nVar), Q2(nVar), D1(nVar)
	REAL :: GRADQ(3, nVar), BGRADQ(nVar), SOURCE(nVAR)
	REAL :: F(nDim,nVar)
	INTEGER:: il, ir, i
	REAL :: x0(3), t

	!PRINT *, "GRMHD Cons2Prim test"
	!PRINT *, "This program tests whether con2prim is the inverse of prim2con"
	
	PRINT *, "GRMHD PDE TEST (nVar=",nVar,"nDim=",nDim,")"

	x0 = (/ 0., 0.1, 0. /)
	t = 0.0
	CALL AlfenWavePrim(x0, t, V0)
	CALL PDEPrim2Cons(Q0, V0)
	
	PRINT *, " Q0 = ", Q0
	
#if 0
	CALL PDECons2Prim(V1, Q0, il)
	PRINT *, "Initial AlfenWave primitives:"
	PRINT *, "V0 =", V0(1:5)
	PRINT *, "AlfenWave Conserved:"
	PRINT *, "Q0  =", Q0(1:5)
	PRINT *, "Reconstructed primitive variables:"
	PRINT *, "V1 =", V1(1:5)
	PRINT *, "Difference:"
	D1 = ABS(V1-V0)
	PRINT *, "ABS(V1-V0) =", D1(1:5)
#endif

	CALL PDEFlux(F(1,:),F(2,:),F(3,:),Q0)
	
	PRINT *, "Fx=",F(1,:)
	PRINT *, "Fy=",F(2,:)
	PRINT *, "Fz=",F(3,:)
	
	CALL DIFFID( x0, t, GRADQ)
	
	PRINT *, "GRADQ(x):", GRADQ(1,:)
	PRINT *, "GRADQ(y):", GRADQ(2,:)
	PRINT *, "GRADQ(z):", GRADQ(3,:)
	
	CALL PDENCP(BGRADQ, Q0, gradQ) 
	SOURCE = -BGRADQ
	
	PRINT *, "SOURCE=", SOURCE
	
	PRINT *, "The FusedSource should be zero, but there is an entry for Sx and tau. Determine why. Examine the derivatives!"
	
	! CALL C MAIN FUNCTION
	! CALL cmain

END ! Program GRMHDStandloneC2P
