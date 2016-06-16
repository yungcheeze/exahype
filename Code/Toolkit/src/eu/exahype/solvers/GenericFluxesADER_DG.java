package eu.exahype.solvers;

import java.io.BufferedWriter;
import java.io.IOException;

public abstract class GenericFluxesADER_DG implements Solver {
  public static final String Identifier = "generic::fluxes::";

  public GenericFluxesADER_DG(int dimensions, int numberOfUnknowns, int numberOfParameters,
      int order) {
    _dimensions = dimensions;
    _numberOfUnknowns = numberOfUnknowns;
    _numberOfParameters = numberOfParameters;
    _order = order;
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
    writer.write(
        "#ifdef EXAHYPE_ENABLE_PROFILER\n  _profiler->start(\"spaceTimePredictor\");\n#endif\n");
    writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
        + "::spaceTimePredictor" + (isLinear() ? "Linear<ncp>" : "Nonlinear<flux>")
        + "( lQi, lFi, lQhi, lFhi, lQhbnd, lFhbnd, luh, dx, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write(
        "#ifdef EXAHYPE_ENABLE_PROFILER\n  _profiler->stop(\"spaceTimePredictor\");\n#endif\n");
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::solutionUpdate(double* luh, const double* const lduh, const double dt) {\n");
    writer
        .write("#ifdef EXAHYPE_ENABLE_PROFILER\n  _profiler->start(\"solutionUpdate\");\n#endif\n");
    writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
        + "::solutionUpdate( luh, lduh, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer
        .write("#ifdef EXAHYPE_ENABLE_PROFILER\n  _profiler->stop(\"solutionUpdate\");\n#endif\n");
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::volumeIntegral(double* lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    writer
        .write("#ifdef EXAHYPE_ENABLE_PROFILER\n  _profiler->start(\"volumeIntegral\");\n#endif\n");
    writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
        + "::volumeIntegral" + (isLinear() ? "Linear" : "Nonlinear")
        + "( lduh, lFhi, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer
        .write("#ifdef EXAHYPE_ENABLE_PROFILER\n  _profiler->stop(\"volumeIntegral\");\n#endif\n");
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::surfaceIntegral(double* lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    writer.write(
        "#ifdef EXAHYPE_ENABLE_PROFILER\n  _profiler->start(\"surfaceIntegral\");\n#endif\n");
    writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
        + "::surfaceIntegral" + (isLinear() ? "Linear" : "Nonlinear")
        + "( lduh, lFhbnd, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer
        .write("#ifdef EXAHYPE_ENABLE_PROFILER\n  _profiler->stop(\"surfaceIntegral\");\n#endif\n");
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::riemannSolver(double* FL, double* FR, const double* const QL, const double* const QR, const double dt, const int normalNonZeroIndex) {\n");
    writer
        .write("#ifdef EXAHYPE_ENABLE_PROFILER\n  _profiler->start(\"riemannSolver\");\n#endif\n");
    writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
        + "::riemannSolver"
        + (isLinear() ? "Linear<eigenvalues, matrixb>" : "Nonlinear<eigenvalues>")
        + "( FL, FR, QL, QR, dt, normalNonZeroIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write("#ifdef EXAHYPE_ENABLE_PROFILER\n  _profiler->stop(\"riemannSolver\");\n#endif\n");
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("double " + projectName + "::" + solverName
        + "::stableTimeStepSize(const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    writer.write(
        "#ifdef EXAHYPE_ENABLE_PROFILER\n  _profiler->start(\"stableTimeStepSize\");\n#endif\n");
    writer.write("   return kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
        + "::stableTimeStepSize<eigenvalues>( luh, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write(
        "#ifdef EXAHYPE_ENABLE_PROFILER\n  _profiler->stop(\"stableTimeStepSize\");\n#endif\n");
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::solutionAdjustment(double *luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {\n");
    writer.write(
        "#ifdef EXAHYPE_ENABLE_PROFILER\n  _profiler->start(\"solutionAdjustment\");\n#endif\n");
    writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
        + "::solutionAdjustment<adjustedSolutionValues>( luh, center, dx, t, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write(
        "#ifdef EXAHYPE_ENABLE_PROFILER\n  _profiler->stop(\"solutionAdjustment\");\n#endif\n");
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::faceUnknownsProlongation(double* lQhbndFine,double* lFhbndFine,const double* lQhbndCoarse,const double* lFhbndCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) {\n");
    writer.write(
        "#ifdef EXAHYPE_ENABLE_PROFILER\n  _profiler->start(\"faceUnknownsProlongation\");\n#endif\n");
    writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
        + "::faceUnknownsProlongation( lQhbndFine, lFhbndFine, lQhbndCoarse, lFhbndCoarse, coarseGridLevel, fineGridLevel, subfaceIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write(
        "#ifdef EXAHYPE_ENABLE_PROFILER\n  _profiler->stop(\"faceUnknownsProlongation\");\n#endif\n");
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::faceUnknownsRestriction(double* lQhbndCoarse,double* lFhbndCoarse,const double* lQhbndFine,const double* lFhbndFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) {\n");
    writer.write(
        "#ifdef EXAHYPE_ENABLE_PROFILER\n  _profiler->start(\"faceUnknownsRestriction\");\n#endif\n");
    writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
        + "::faceUnknownsRestriction( lQhbndCoarse, lFhbndCoarse, lQhbndFine, lFhbndFine, coarseGridLevel, fineGridLevel, subfaceIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write(
        "#ifdef EXAHYPE_ENABLE_PROFILER\n  _profiler->stop(\"faceUnknownsRestriction\");\n#endif\n");
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::volumeUnknownsProlongation(double* luhFine, const double* luhCoarse, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {\n");
    writer.write(
        "#ifdef EXAHYPE_ENABLE_PROFILER\n  _profiler->start(\"volumeUnknownsProlongation\");\n#endif\n");
    writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
        + "::volumeUnknownsProlongation( luhFine, luhCoarse, coarseGridLevel, fineGridLevel, subcellIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write(
        "#ifdef EXAHYPE_ENABLE_PROFILER\n  _profiler->stop(\"volumeUnknownsProlongation\");\n#endif\n");
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::volumeUnknownsRestriction(double* luhCoarse, const double* luhFine, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {\n");
    writer.write(
        "#ifdef EXAHYPE_ENABLE_PROFILER\n  _profiler->start(\"volumeUnknownsRestriction\");\n#endif\n");
    writer.write("  kernels::aderdg::generic::" + (isFortran() ? "fortran" : "c")
        + "::volumeUnknownsRestriction( luhCoarse, luhFine, coarseGridLevel, fineGridLevel, subcellIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write(
        "#ifdef EXAHYPE_ENABLE_PROFILER\n  _profiler->stop(\"volumeUnknownsRestriction\");\n#endif\n");
    writer.write("}\n");
    writer.write("\n\n\n");
  }

  protected int _dimensions;
  protected int _numberOfUnknowns;
  protected int _numberOfParameters;
  protected int _order;
}
