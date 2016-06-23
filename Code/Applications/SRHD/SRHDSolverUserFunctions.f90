 ! 
 ! This file is part of the ExaHyPE project.
 ! Copyright (c) 2016  http://exahype.eu
 ! All rights reserved.
 !
 ! The project has received funding from the European Union's Horizon 
 ! 2020 research and innovation programme under grant agreement
 ! No 671698. For copyrights and licensing, please consult the webpage.
 !
 ! Released under the BSD 3 Open Source License.
 ! For the full license text, see LICENSE.txt
 !  
 
 SUBROUTINE MinimumTreeDepth(depth)
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE 
    ! Argument list 
    INTEGER, INTENT(OUT)              :: depth        ! maximal depth of tree recursion
    
    depth = 4
    
END SUBROUTINE MinimumTreeDepth 

SUBROUTINE HasToAdjustSolution(time, refine)
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE 
    ! Argument list 
    REAL   , INTENT(IN)               :: time        ! 

    LOGICAL, INTENT(OUT)              :: refine      ! 
    
    IF(time<0.000000001) THEN
      refine = .TRUE.
    ELSE
      refine = .FALSE.
    ENDiF
    
END SUBROUTINE HasToAdjustSolution


SUBROUTINE AdjustedSolutionValues(x, w, t, dt, Q)
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)               :: x(2)        ! 
    REAL, INTENT(IN)               :: w           ! 
    REAL, INTENT(IN)               :: t           ! 
    REAL, INTENT(IN)               :: dt          ! 

    REAL, INTENT(OUT)              :: Q(5)        ! 
    
    REAL :: GAMMA

    !PRINT *, ' ---------------------------------------- ' 
    !PRINT *, x
    !PRINT *, ' ---------------------------------------- ' 
    !PRINT *, w
    !PRINT *, ' ---------------------------------------- ' 
    !PRINT *, t
    !PRINT *, ' ---------------------------------------- ' 
    !PRINT *, dt
    !PRINT *, ' ---------------------------------------- ' 
    !PRINT *, Q
    !PRINT *, ' ---------------------------------------- ' 
    
    GAMMA = 1.4
    Q(1) = 1.
    Q(2) = 0.
    Q(3) = 0.
    Q(4) = 0.
    !Q(5) = EXP((x(1)-0.5)*(x(1)-0.5) + (x(2)-0.5)*(x(2)-0.5))
    Q(5) = 1. / (GAMMA - 1) + exp(-((x(1) - 0.5) * (x(1) - 0.5) + (x(2) - 0.5) * (x(2) - 0.5)) / (0.05 * 0.05)) * 1.0e-3

    
END SUBROUTINE AdjustedSolutionValues


SUBROUTINE Flux(Q, F)
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)              :: Q(5)        ! 

    REAL, INTENT(OUT)             :: F(5, 2)        ! 

    REAL :: p, irho  

    !PRINT *, ' ---------------------------------------- ' 
    !PRINT *, F
    !PRINT *, ' ---------------------------------------- ' 
    !PRINT *, Q
    !PRINT *, ' ---------------------------------------- ' 
    
    irho = 1.0/Q(1)
    !p = (EQNgamma-1)*( Q(5) - 0.5*SUM(Q(2:4)**2)*irho )
    p  = (0.4       )*( Q(5) - 0.5*SUM(Q(2:4)**2)*irho )
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
    
   
END SUBROUTINE Flux


SUBROUTINE Eigenvalues(Q, nnzi, lambda)
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE 
    ! Argument list 
    REAL, INTENT(IN)              :: Q(5)        ! 
    INTEGER, INTENT(IN)           :: nnzi        ! 

    REAL, INTENT(OUT)             :: lambda(5)   ! 

    REAL :: p, u, c, nv(2)
    
    nv(:) = 0
    nv(nnzi) = 1
    
    !
    u = ( Q(2)*nv(1) + Q(3)*nv(2) + Q(4)*nv(3) )/Q(1)       ! normal velocity 
    p  = (0.4       )*( Q(5) - 0.5*SUM(Q(2:4)**2)/Q(1) )    ! fluid pressure 
    c  = SQRT(1.4     *p/Q(1))                              ! sound speed
    !
    Lambda = (/ u-c, u, u, u, u+c /)                        ! The eigenvalues of the Euler equations 
   
   
END SUBROUTINE Eigenvalues



