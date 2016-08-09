package eu.exahype.solvers;

public final class GenericFluxesNonlinearADER_DGinFortran extends GenericFluxesADER_DG {
  public static final String Identifier = GenericFluxesNonlinearADER_DGinC.Identifier;

  public GenericFluxesNonlinearADER_DGinFortran(int dimensions, int numberOfUnknowns,
      int numberOfParameters, int order, boolean enableProfiler, boolean hasConstants) {
    super(dimensions, numberOfUnknowns, numberOfParameters, order, enableProfiler, hasConstants);
  }

  @Override
  public final boolean isLinear() {
    return false;
  }

  @Override
  public final boolean isFortran() {
    return true;
  }

  public void writeUserImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException {
    Helpers.writeMinimalADERDGSolverUserImplementation(solverName, writer, projectName,
        _numberOfUnknowns, _numberOfParameters, _order, _hasConstants);

    int digits = String.valueOf(_numberOfUnknowns + _numberOfParameters).length();

    writer.write("bool " + projectName + "::" + solverName
        + "::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t) {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("  return false;\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {\n");
    writer.write("  // Dimensions             = " + _dimensions + "\n");
    writer.write(
        "  // Number of variables    = " + Integer.toString(_numberOfUnknowns + _numberOfParameters)
            + " (#unknowns + #parameters)\n");
    writer.write("  // @todo Please implement\n");
    for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
      writer.write("  Q[" + String.format("%" + digits + "d", i) + "] = 0.0;\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("exahype::solvers::Solver::RefinementControl " + projectName + "::" + solverName
        + "::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("  return exahype::solvers::Solver::RefinementControl::Keep;\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("\n\n\n");
    writer.write("//************************************************* \n");
    writer.write("//for FORTRAN kernels the fluxes and eigenvalues \n");
    writer.write("//have to be implemented in the file ./PDE.f90 \n");
    writer.write("//************************************************* \n");
    writer.write("\n\n\n");
  }

  public void writeUserPDE(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    // @todo Implement
    int digits = String.valueOf(_numberOfUnknowns + _numberOfParameters).length();

    writer.write("SUBROUTINE PDEEigenvalues(Lambda,Q,nv) \n");
    writer.write("  USE typesDef, ONLY : nVar, d \n");
    writer.write("  USE, INTRINSIC :: ISO_C_BINDING \n");
    writer.write("  IMPLICIT NONE \n");
    writer.write("  ! Argument list  \n");
    writer.write("  REAL, INTENT(IN)  :: Q(nVar), nv(d)  \n");
    writer.write("  REAL, INTENT(OUT) :: Lambda(nVar)  \n");
    writer.write("  ! Local variables  \n");
    writer.write("  !\n");
    writer.write("  !@todo Please implement\n");
    for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
      writer.write("  Lambda(" + String.format("%" + digits + "d", i + 1) + ") = 0.0\n");
    }
    writer.write("  !\n");
    writer.write("END SUBROUTINE PDEEigenvalues\n");

    writer.write(" \n\n\n");

    writer.write("SUBROUTINE PDEFlux(F,Q) \n");
    writer.write("  USE typesDef, ONLY : nVar, d \n");
    writer.write("  USE, INTRINSIC :: ISO_C_BINDING \n");
    writer.write("  IMPLICIT NONE \n");
    writer.write("  ! Argument list  \n");
    writer.write("  REAL, INTENT(IN)  :: Q(nVar) \n");
    writer.write("  REAL, INTENT(OUT) :: F(nVar,d) \n");
    writer.write("  ! Local variables  \n");
    writer.write("  !\n");
    writer.write("  !@todo Please implement\n");
    writer.write("  !\n");
    for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
      writer.write("  F(" + String.format("%" + digits + "d", i + 1) + ", 1) = 0.0\n");
    }
    writer.write("  !\n");
    for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
      writer.write("  F(" + String.format("%" + digits + "d", i + 1) + ", 2) = 0.0\n");
    }
    if (_dimensions == 3) {
      writer.write("  !\n");
      for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
        writer.write("  F(" + String.format("%" + digits + "d", i + 1) + ", 3) = 0.0\n");
      }
    }
    writer.write("  !\n");
    writer.write("END SUBROUTINE PDEFlux \n");

    writer.write(" \n\n\n");

    writer.write("SUBROUTINE PDENCP(BgradQ,Q,gradQ) \n");
    writer.write("  ! not used  \n");
    writer.write("  PRINT *, 'PDENCP is not used for nonlinear-solvers' \n");
    writer.write("  CALL EXIT  \n");
    writer.write("  !\n");
    writer.write("END SUBROUTINE PDENCP \n");

    writer.write(" \n\n\n");

    writer.write("SUBROUTINE PDEMatrixB(Bn,Q,nv) \n");
    writer.write("  ! not used  \n");
    writer.write("  PRINT *, 'PDEMatrixB is not used for nonlinear-solvers' \n");
    writer.write("  CALL EXIT  \n");
    writer.write("  !\n");
    writer.write("END SUBROUTINE PDEMatrixB \n");
  }

  public void writeTypesDef(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    // @todo Implement

    writer.write("  MODULE typesDef \n");
    writer.write("    IMPLICIT NONE  \n");
    writer.write("    PUBLIC  \n");
    writer.write("    ! \n");
    writer.write(
        "    ! ================================== This part of the typesDef can be modified by the user.  ==================================  \n");
    writer.write("    ! \n");
    writer.write(
        "    INTEGER, PARAMETER             :: d = 3                               ! This is the maximum number of space dimensions we want to deal with in our heads. !! NEVER change this parameter, unless you are bold and want to solve the Boltzmann equation !!  \n");
    writer.write("    INTEGER, PARAMETER             :: N = " + _order
        + "                               ! Polynomial degree of our approximation in space and time  \n");
    writer.write("    INTEGER, PARAMETER             :: nDim = " + _dimensions
        + "                            ! The number of space dimensions that we actually want to simulate  \n");
    writer.write("    INTEGER, PARAMETER             :: nVar = " + _numberOfUnknowns
        + "                            ! The number of variables of the PDE system  \n");
    writer.write("    INTEGER, PARAMETER             :: nDOF(0:3) = (/ " + (_order + 1) + ", "
        + (_order + 1) + ", " + (_order + 1) + ", " + (_order + 1)
        + " /)                           ! The number of degrees of freedom in space and time  \n");
    writer.write("     \n");
    writer.write(
        "    DOUBLE PRECISION, PARAMETER    :: wGPN(N+1)     = (/ 0.1739274225687273, 0.3260725774312732, 0.3260725774312732, 0.1739274225687273 /) \n");
    writer.write(
        "    DOUBLE PRECISION, PARAMETER    :: xiGPN(N+1)    = (/ 0.0694318442029737, 0.3300094782075719, 0.6699905217924281, 0.9305681557970262 /) \n");
    writer.write(
        "    DOUBLE PRECISION, PARAMETER    :: F0(N+1)       = (/ 1.526788125457e+00, -8.136324494869e-01, 4.007615203117e-01, -1.139171962820e-01 /) \n");
    writer.write(
        "    DOUBLE PRECISION, PARAMETER    :: FLcoeff(N+1)  = (/ 1.52678812545727, -0.813632449486927, 0.400761520311650, -0.113917196281990 /) \n");
    writer.write(
        "    DOUBLE PRECISION, PARAMETER    :: FRcoeff(N+1)  = (/ -0.113917196281990, 0.400761520311651, -0.813632449486928, 1.52678812545727 /) \n");
    writer.write(
        "    DOUBLE PRECISION, PARAMETER    :: dudx(N+1,N+1) = reshape( (/ -6.66400047270456, -1.51511522959847, 0.657396448516548, -1.16125633832453, 9.72030883137039, -0.768828784446417, -2.94134046256143, 4.21756469699036, -4.21756469699036, 2.94134046256143, 0.768828784446416, -9.72030883137039, 1.16125633832453, -0.657396448516549, 1.51511522959847, 6.66400047270456 /), (/N+1,N+1 /) ) \n");
    writer.write(
        "    DOUBLE PRECISION, PARAMETER    :: Kxi(N+1,N+1)  = reshape( (/ -1.15905242621428, 1.69062826161229, -0.733550157264387, 0.201974321866383, -0.494037528020548, -0.250693983347796, 0.959090465730300, -0.214358954361956, 0.214358954361956, -0.959090465730300, 0.250693983347796, 0.494037528020548, -0.201974321866383, 0.733550157264387, -1.69062826161229, 1.15905242621428 /), (/N+1,N+1 /) ) \n");
    writer.write(
        "    DOUBLE PRECISION, PARAMETER    :: iK1(N+1,N+1)  = reshape( (/ 0.546435362419645, 1.01885331677130, 1.02401050669309, 0.974005058264396, -0.144326183293257, 0.584759972857323, 1.00074377855320, 1.02401050669309, 0.101462359828863, -0.170263724267844, 0.584759972857323, 1.01885331677130, -6.687578310368468E-002, 0.101462359828862, -0.144326183293257, 0.546435362419644 /), (/N+1,N+1 /) ) \n");
    writer.write("     \n");
    writer.write("   \n");
    writer.write("    TYPE tFace \n");
    writer.write(
        "      DOUBLE PRECISION, POINTER    :: qL(:,:,:), qR(:,:,:)                ! pointer to left and right boundary-extrapolated state vector  \n");
    writer.write(
        "      DOUBLE PRECISION, POINTER    :: FL(:,:,:), FR(:,:,:)                ! pointer to left and right boundary-extrapolated flux vector  \n");
    writer.write(
        "      INTEGER          :: Left, Right                         ! pointer to left and right element  \n");
    writer.write(
        "      DOUBLE PRECISION             :: nv(d)                               ! face normal vector  \n");
    writer.write("    END TYPE       \n");
    writer.write("    TYPE(tFace), POINTER :: Face(:)  \n");
    writer.write("  END MODULE typesDef  \n");
  }
}
