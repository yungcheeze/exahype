! GRMHD ID Copied

RECURSIVE SUBROUTINE AlfenWavePrim(x, t, Q)
    ! Computes the AlfenWave conserved variables (Q) at a given time t.
    ! Use it ie. with t=0 for initial data
    ! Use it for any other time ie. for comparison
    
    ! GRID FOR ALFENWAVE:
    !     dimension const                = 2
    !     width                          = 1.0, 0.3
    !     offset                         = 0.0, 0.0
    !     end-time                       = 2.1
    !
    !  maximum-mesh-size              = 0.04

    USE, INTRINSIC :: ISO_C_BINDING
    USE Parameters, ONLY : nVar, nDim, gamma
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: t
    REAL, INTENT(IN)               :: x(nDim)        ! 
    REAL, INTENT(OUT)              :: Q(nVar)        ! 
    REAL, PARAMETER                :: t_offset = 1.0 ! time offset

    REAL :: rho0, p0, eta, B0, hh, tempaa, tempab, tempac, va2, vax
    REAL :: V(nVar), BV(3), VV(3), Pi = ACOS(-1.0)

    rho0 = 1.
    p0   = 1.
    eta  = 1.
    B0   = 1.0 
    !
    hh = 1.0 + gamma / ( gamma - 1.0) * p0 / rho0
    tempaa = rho0 * hh + B0**2 * ( 1.0 + eta**2)
    tempab = 2.0 * eta * B0**2 / tempaa
    tempac = 0.5 * ( 1.0 + sqrt ( 1.0 - tempab**2))
    va2 = b0**2 / ( tempaa * tempac)
    vax = sqrt ( va2)
    !
    BV(1) = B0
    BV(2) = eta * B0 * COS(2*Pi*( x(1) - vax*(t-t_offset)))
    BV(3) = eta * B0 * SIN(2*Pi*( x(1) - vax*(t-t_offset)))
    !
    VV(1)   = 0.0
    VV(2:3) = - vax * BV(2:3) / B0
    !
    ! Now convert to conservative variables
    !
    V(1:9) = (/ rho0, VV(1:3), p0, BV(1:3), 0.0 /) ! psi
    !            lapse, shift(3),    gamma11, ...  gamma22, -, gamma33
    V(10:19) = (/ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0 /)
    
    
    !CALL PDEPrim2Cons(Q,V)
    ! JUST COPY PRIMITIVES TO CONSERVED
    Q = V
END SUBROUTINE AlfenWavePrim

