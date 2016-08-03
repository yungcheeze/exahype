package eu.exahype.solvers;

import java.io.BufferedWriter;
import java.io.IOException;

public abstract class GenericFluxesADER_DG implements Solver {
  public static final String Identifier = "generic::fluxes::";

  public GenericFluxesADER_DG(int dimensions, int numberOfUnknowns, int numberOfParameters,
      int order, boolean enableProfiler) {
    _dimensions = dimensions;
    _numberOfUnknowns = numberOfUnknowns;
    _numberOfParameters = numberOfParameters;
    _order = order;
    _enableProfiler = enableProfiler;
  }

  public abstract boolean isLinear();

  public abstract boolean isFortran();

  @Override
  public final void writeHeader(BufferedWriter writer, String solverName, String projectName)
      throws IOException {
    Helpers.writeMinimalADERDGSolverHeader(solverName, writer, projectName);

    writer.write("  private:\n");
    writer.write(
        "    static void eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda);\n");
    writer.write("    static void flux(const double* const Q, double** F);\n");
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
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::spaceTimePredictor( double* lQi, double* lFi, double* lQhi, double* lFhi, double* lQhbnd, double* lFhbnd, const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx, const double dt ) {\n");
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"spaceTimePredictor\");\n");
    }
    writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
        + "::spaceTimePredictor" + (isLinear() ? "Linear<ncp>" : "Nonlinear<flux>")
        + "( lQi, lFi, lQhi, lFhi, lQhbnd, lFhbnd, luh, dx, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
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
        + "::riemannSolver(double* FL, double* FR, const double* const QL, const double* const QR, const double dt, const int normalNonZeroIndex) {\n");
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"riemannSolver\");\n");
    }
    writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
        + "::riemannSolver"
        + (isLinear() ? "Linear<eigenvalues, matrixb>" : "Nonlinear<eigenvalues>")
        + "( FL, FR, QL, QR, dt, normalNonZeroIndex, getNumberOfVariables(), getNumberOfParameters(), getNodesPerCoordinateAxis() );\n");
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
        + "::stableTimeStepSize(const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"stableTimeStepSize\");\n");
    }
    writer.write("  double d = kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
        + "::stableTimeStepSize<eigenvalues>( luh, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
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
