MODULE typesDef
  IMPLICIT NONE 
  PUBLIC 
  !
  ! ================================== This part of the typesDef can be modified by the user.  ================================== 
  !
  INTEGER, PARAMETER             :: d = 3                               ! This is the maximum number of space dimensions we want to deal with in our heads. !! NEVER change this parameter, unless you are bold and want to solve the Boltzmann equation !! 
  INTEGER, PARAMETER             :: N = 3                               ! Polynomial degree of our approximation in space and time 
  INTEGER, PARAMETER             :: nDim = 3                            ! The number of space dimensions that we actually want to simulate 
  DOUBLE PRECISION, PARAMETER    :: CFL = 0.9                           ! The Courant-Friedrichs-Lewy number < 1 
  INTEGER, PARAMETER             :: nVar = 5                            ! The number of variables of the PDE system 
  INTEGER, PARAMETER             :: nDOF(0:3) = (/ 4, 4, 4, 4 /)                           ! The number of degrees of freedom in space and time 
  DOUBLE PRECISION, PARAMETER    :: EQNgamma = 1.4                           ! The Courant-Friedrichs-Lewy number < 1 
  
  DOUBLE PRECISION, PARAMETER    :: wGPN(N+1)     = (/ 0.1739274225687273, 0.3260725774312732, 0.3260725774312732, 0.1739274225687273 /)
  DOUBLE PRECISION, PARAMETER    :: xiGPN(N+1)    = (/ 0.0694318442029737, 0.3300094782075719, 0.6699905217924281, 0.9305681557970262 /)
  DOUBLE PRECISION, PARAMETER    :: F0(N+1)       = (/ 1.526788125457e+00, -8.136324494869e-01, 4.007615203117e-01, -1.139171962820e-01 /)
  DOUBLE PRECISION, PARAMETER    :: FLcoeff(N+1)  = (/ 1.52678812545727, -0.813632449486927, 0.400761520311650, -0.113917196281990 /)
  DOUBLE PRECISION, PARAMETER    :: FRcoeff(N+1)  = (/ -0.113917196281990, 0.400761520311651, -0.813632449486928, 1.52678812545727 /)
  DOUBLE PRECISION, PARAMETER    :: Kxi(N+1,N+1)  = (/ -1.15905242621428, 1.69062826161229, -0.733550157264387, 0.201974321866383, -0.494037528020548, -0.250693983347796, 0.959090465730300, -0.214358954361956, 0.214358954361956, -0.959090465730300, 0.250693983347796, 0.494037528020548, -0.201974321866383, 0.733550157264387, -1.69062826161229, 1.15905242621428 /)
  DOUBLE PRECISION, PARAMETER    :: iK1(N+1,N+1)  = (/ 0.546435362419645, 1.01885331677130, 1.02401050669309, 0.974005058264396, -0.144326183293257, 0.584759972857323, 1.00074377855320, 1.02401050669309, 0.101462359828863, -0.170263724267844, 0.584759972857323, 1.01885331677130, -6.687578310368468E-002, 0.101462359828862, -0.144326183293257, 0.546435362419644 /)
  
  !
  ! ==================================           Do NOT change the stuff below                 ==================================
  !
  ! The following variables contain important information about the numerical method. Do NOT change.  
  !
  !DOUBLE PRECISION, PARAMETER    :: PNPMTable(0:9) = (/ 1.0, 0.33, 0.17, 0.1, 0.069, 0.045,  0.038, 0.03, 0.02, 0.015 /)   ! maximum CFL numbers for PNPM schemes according to von Neumann stability analysis and experience     
  !DOUBLE PRECISION               :: xin(N+1)                            ! The nodes used for our basis (can in principle be different from the Gauss-Legendre nodes, like for example also the bad Gauss-Lobatto nodes) 
  !DOUBLE PRECISION               :: MM(N+1,N+1), iMM(N+1,N+1)           ! Element mass-matrix and its inverse 
  !DOUBLE PRECISION                :: dudx(N+1,N+1)         ! Element stiffness matrix and discrete derivative operator 
  !DOUBLE PRECISION               :: FLm(N+1,N+1), FLp(N+1,N+1)          ! Left flux matrices 
  !DOUBLE PRECISION               :: FRm(N+1,N+1), FRp(N+1,N+1)          ! Right flux matrices 
  !DOUBLE PRECISION               :: FLcoeff(N+1), FRcoeff(N+1)          
  !DOUBLE PRECISION               :: F1(N+1,N+1)                ! Time flux matrices 
  !DOUBLE PRECISION               :: K1(N+1,N+1)           ! F1 - Ktau 
  !INTEGER            :: dn(d)                               ! number of direct neighbors in each dimension 
  
  ! Stuff related to the problem setup, mesh and output 
  !INTEGER            :: IMAX, JMAX, KMAX, NMAX              ! The number of cells in each space dimension & max. number of time steps 
  !INTEGER            :: nElem, nFace, nNode                 ! The number of elements, faces and nodes 
  !INTEGER            :: timestep                            ! the number of the current time step 
  !DOUBLE PRECISION               :: xL(d), xR(d)                        ! computational domain 
  !DOUBLE PRECISION               :: dx(d), dt                           ! The vector of mesh spacings in each dimension and the time step   
  !DOUBLE PRECISION               :: time, tend                          ! current time  and final time 
  !DOUBLE PRECISION               :: tio, dtio                           ! output time and output time interval 
  !DOUBLE PRECISION, POINTER      :: x(:,:)                              ! the node coordinates (nDim, nNode) 
  !INTEGER, PARAMETER :: nVtx = 2**nDim, nFac = 2*nDim       ! number of vertices and faces per element 
  !INTEGER, POINTER   :: tri(:,:)                            ! connectivity from the element to the nodes 
  !INTEGER, POINTER   :: Element2Face(:,:)                   ! connectivity from each element to its faces 
  !CHARACTER(LEN=200) :: BaseFile                            ! Basic filename to write the results 
  !DOUBLE PRECISION               :: SubOutputMatrix((N+1)**d,(N+1)**d)  ! Matrix needed for the plotting of the results on a fine subgrid 
  !INTEGER            :: subtri(2**d,N**d)                   ! subcell connectivity (for fine output) 
  !DOUBLE PRECISION               :: allsubxi(d,(N+1)**d)                ! subnodes (for fine output) 
  
  ! Some diagnostics                                        ! 
  !DOUBLE PRECISION               :: tCPU1, tCPU2                        ! CPU times 
  !INTEGER(8)         :: TEU                                 ! total element updates 
  !
  
  TYPE tFace
    DOUBLE PRECISION, POINTER    :: qL(:,:,:), qR(:,:,:)                ! pointer to left and right boundary-extrapolated state vector 
    DOUBLE PRECISION, POINTER    :: FL(:,:,:), FR(:,:,:)                ! pointer to left and right boundary-extrapolated flux vector 
    INTEGER          :: Left, Right                         ! pointer to left and right element 
    DOUBLE PRECISION             :: nv(d)                               ! face normal vector 
  END TYPE      
  TYPE(tFace), POINTER :: Face(:) 
  !
  ! The main variables of the ADER-DG scheme 
  !
  ! Important info and parameters concerning the governing PDE system 
  !
  
  !
END MODULE typesDef 
    
    
    
    
