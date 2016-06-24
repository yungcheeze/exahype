  MODULE typesDef 
    IMPLICIT NONE  
    PUBLIC  
    ! 
    ! ================================== This part of the typesDef can be modified by the user.  ==================================  
    ! 
    INTEGER, PARAMETER             :: d = 3                               ! This is the maximum number of space dimensions we want to deal with in our heads. !! NEVER change this parameter, unless you are bold and want to solve the Boltzmann equation !!  
    INTEGER, PARAMETER             :: N = 3                               ! Polynomial degree of our approximation in space and time  
    INTEGER, PARAMETER             :: nDim = 3                            ! The number of space dimensions that we actually want to simulate  
    INTEGER, PARAMETER             :: nVar = 5                            ! The number of variables of the PDE system  
    INTEGER, PARAMETER             :: nDOF(0:3) = (/ 4, 4, 4, 4 /)                           ! The number of degrees of freedom in space and time
  END MODULE typesDef  
