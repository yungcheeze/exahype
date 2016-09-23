package eu.exahype.solvers;

import java.io.BufferedWriter;
import java.io.IOException;

public abstract class GenericFluxesADER_DG implements Solver {
  public static final String Identifier = "generic::fluxes::";
  protected boolean _hasConstants;

  public GenericFluxesADER_DG(int dimensions, int numberOfUnknowns, int numberOfParameters,
      int order, boolean enableProfiler, boolean hasConstants) {
    _dimensions = dimensions;
    _numberOfUnknowns = numberOfUnknowns;
    _numberOfParameters = numberOfParameters;
    _order = order;
    _enableProfiler = enableProfiler;
    _hasConstants   = hasConstants;
  }

  public abstract boolean isLinear();

  public abstract boolean isFortran();

  @Override
  public final void writeHeader(BufferedWriter writer, String solverName, String projectName)
      throws IOException {
    Helpers.writeMinimalADERDGSolverHeader(solverName, writer, projectName, _hasConstants );

    writer.write("  private:\n");
    writer.write("    void init();\n");
    writer.write("    static void eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda);\n");
    writer.write("    static void flux(const double* const Q, double** F);\n");
    writer.write("    static void source(const double* const Q, double* S);\n");
    writer.write("    static void boundaryValues(const double* const x,const double t, const int faceIndex, const int normalNonZero, const double * const fluxIn, const double* const stateIn, double *fluxOut, double* stateOut);\n");
    writer.write(
        "    static void adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q);\n");

    if (isLinear()) {
      writer.write(
          "    static void ncp(const double* const Q, const double* const gradQ, double* BgradQ);\n");
      writer.write(
          "    static void matrixb(const double* const Q, const int normalNonZero, double* Bn);\n");
    }

    writer.write("};\n\n\n");
  }

  @Override
  public final void writeGeneratedImplementation(BufferedWriter writer, String solverName,
      String projectName) throws IOException {

    writer.write("// ==============================================\n");
    writer.write("// Please do not change the implementations below\n");
    writer.write("// =============================---==============\n");
    writer.write("#include \"" + solverName + ".h\"\n");
    writer.write("#include \"kernels/aderdg/generic/Kernels.h\"\n");
    writer.write("\n\n");

    // constructor
    if (_hasConstants) {
      writer.write(projectName + "::" + solverName + "::" + solverName + "(double maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping, std::unique_ptr<exahype::profilers::Profiler> profiler, exahype::Parser::ParserView constants):\n");
    }
    else {
      writer.write(projectName + "::" + solverName + "::" + solverName + "(double maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping, std::unique_ptr<exahype::profilers::Profiler> profiler):\n");
    }


    writer.write("  exahype::solvers::ADERDGSolver("
        + "\""+solverName+"\", "+_numberOfUnknowns+" /* numberOfUnknowns */, "
        +_numberOfParameters+"/* numberOfParameters */, " + (_order + 1) +
        " /* nodesPerCoordinateAxis */, maximumMeshSize, timeStepping, " +
        "std::move(profiler)) {\n");
    writer.write("  init();\n");
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::spaceTimePredictor(double* lQhbnd,double* lFhbnd,double** tempSpaceTimeUnknowns,double** tempSpaceTimeFluxUnknowns,double*  tempUnknowns,double*  tempFluxUnknowns,const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& dx,const double dt) {\n");
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"spaceTimePredictor\");\n");
    }
    if (isLinear()) {
      writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
          + "::spaceTimePredictorLinear<ncp>"
          + "( lQi, lFi, lQhi, lFhi, lQhbnd, lFhbnd, luh, dx, dt, getNumberOfVariables(), getNumberOfParameters(), getNodesPerCoordinateAxis() );\n");
    } else {
      writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
          + "::spaceTimePredictorNonlinear<flux, source>"
          + "( lQhbnd, lFhbnd, tempSpaceTimeUnknowns, tempSpaceTimeFluxUnknowns, tempUnknowns, tempFluxUnknowns, luh, dx, dt, getNumberOfVariables(), getNumberOfParameters(), getNodesPerCoordinateAxis() );\n");
    }
    if (_enableProfiler) {
      writer.write("  _profiler->stop(\"spaceTimePredictor\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::solutionUpdate(double* luh, const double* const lduh, const double dt) {\n");
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"solutionUpdate\");\n");
    }
    writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
        + "::solutionUpdate( luh, lduh, dt, getNumberOfVariables(), getNumberOfParameters(), getNodesPerCoordinateAxis() );\n");
    if (_enableProfiler) {
      writer.write("  _profiler->stop(\"solutionUpdate\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::volumeIntegral(double* lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"volumeIntegral\");\n");
    }
    writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
        + "::volumeIntegral" + (isLinear() ? "Linear" : "Nonlinear")
        + "( lduh, lFhi, dx, getNumberOfVariables(), getNumberOfParameters(), getNodesPerCoordinateAxis() );\n");
    if (_enableProfiler) {
      writer.write("  _profiler->stop(\"volumeIntegral\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::surfaceIntegral(double* lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"surfaceIntegral\");\n");
    }
    writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
        + "::surfaceIntegral" + (isLinear() ? "Linear" : "Nonlinear")
        + "( lduh, lFhbnd, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    if (_enableProfiler) {
      writer.write("  _profiler->stop(\"surfaceIntegral\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::riemannSolver(double* FL, double* FR, const double* const QL, const double* const QR, double* tempFaceUnknownsArray, double** tempStateSizedVectors, double** tempStateSizedSquareMatrices, const double dt, const int normalNonZeroIndex) {\n");
    
    writer.write("  assertion2(normalNonZeroIndex>=0,dt,normalNonZeroIndex);\n");
    writer.write("  assertion2(normalNonZeroIndex<DIMENSIONS,dt,normalNonZeroIndex);\n");
    
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"riemannSolver\");\n");
    }
    if (isLinear()) {
        writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
                + "::riemannSolverLinear<eigenvalues, matrixb>"
                + "( FL, FR, QL, QR, dt, normalNonZeroIndex, getNumberOfVariables(), getNumberOfParameters(), getNodesPerCoordinateAxis() );\n");
    } else {
        writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
                + "::riemannSolverNonlinear<eigenvalues>"
                + "( FL, FR, QL, QR, tempFaceUnknownsArray, tempStateSizedVectors, tempStateSizedSquareMatrices, dt, normalNonZeroIndex, getNumberOfVariables(), getNumberOfParameters(), getNodesPerCoordinateAxis() );\n");
    }
    if (_enableProfiler) {
      writer.write("  _profiler->stop(\"riemannSolver\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");

    // boundaryConditions
    writer.write("void " + projectName + "::" + solverName
            + "::boundaryConditions(double* fluxOut,double* stateOut,const double* const fluxIn,const double* const stateIn,const tarch::la::Vector<DIMENSIONS, double>& cellCentre,const tarch::la::Vector<DIMENSIONS,double>& cellSize,const double t,const double dt,const int faceIndex,const int normalNonZero) {\n");
    if (_enableProfiler) {
        writer.write("  _profiler->start(\"boundaryConditions\");\n");
    }
    writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
            + "::boundaryConditions<boundaryValues>"
            + "( fluxOut, stateOut, fluxIn, stateIn, cellCentre, cellSize, t, dt, faceIndex, normalNonZero, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    if (_enableProfiler) {
        writer.write("  _profiler->stop(\"boundaryConditions\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
    
    writer.write("double " + projectName + "::" + solverName
        + "::stableTimeStepSize(const double* const luh,double* tempEigenvalues,const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"stableTimeStepSize\");\n");
    }
    writer.write("  double d = kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
        + "::stableTimeStepSize<eigenvalues>( luh, tempEigenvalues, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    if (_enableProfiler) {
      writer.write("  _profiler->stop(\"stableTimeStepSize\");\n");
    }
    writer.write("  return d;\n");
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::solutionAdjustment(double *luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {\n");
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"solutionAdjustment\");\n");
    }
    writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
        + "::solutionAdjustment<adjustedSolutionValues>( luh, center, dx, t, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    if (_enableProfiler) {
      writer.write("  _profiler->stop(\"solutionAdjustment\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::faceUnknownsProlongation(double* lQhbndFine,double* lFhbndFine,const double* lQhbndCoarse,const double* lFhbndCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) {\n");
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"faceUnknownsProlongation\");\n");
    }
    writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
        + "::faceUnknownsProlongation( lQhbndFine, lFhbndFine, lQhbndCoarse, lFhbndCoarse, coarseGridLevel, fineGridLevel, subfaceIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    if (_enableProfiler) {
      writer.write("  _profiler->stop(\"faceUnknownsProlongation\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::faceUnknownsRestriction(double* lQhbndCoarse,double* lFhbndCoarse,const double* lQhbndFine,const double* lFhbndFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) {\n");
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"faceUnknownsRestriction\");\n");
    }
    writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
        + "::faceUnknownsRestriction( lQhbndCoarse, lFhbndCoarse, lQhbndFine, lFhbndFine, coarseGridLevel, fineGridLevel, subfaceIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    if (_enableProfiler) {
      writer.write("  _profiler->stop(\"faceUnknownsRestriction\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::volumeUnknownsProlongation(double* luhFine, const double* luhCoarse, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {\n");
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"volumeUnknownsProlongation\");\n");
    }
    writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
        + "::volumeUnknownsProlongation( luhFine, luhCoarse, coarseGridLevel, fineGridLevel, subcellIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    if (_enableProfiler) {
      writer.write("  _profiler->stop(\"volumeUnknownsProlongation\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::volumeUnknownsRestriction(double* luhCoarse, const double* luhFine, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {\n");
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"volumeUnknownsRestriction\");\n");
    }
    writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
        + "::volumeUnknownsRestriction( luhCoarse, luhFine, coarseGridLevel, fineGridLevel, subcellIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    if (_enableProfiler) {
      writer.write("  _profiler->stop(\"volumeUnknownsRestriction\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
  }

  protected int _dimensions;
  protected int _numberOfUnknowns;
  protected int _numberOfParameters;
  protected int _order;
  protected boolean _enableProfiler;
}
