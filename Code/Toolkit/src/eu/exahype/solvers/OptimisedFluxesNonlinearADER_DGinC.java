package eu.exahype.solvers;

public class OptimisedFluxesNonlinearADER_DGinC implements Solver {
  public static final String Identifier = "optimised::fluxes::nonlinear";

  private int _dimensions;
  private int _numberOfUnknowns;
  private int _numberOfParameters;
  private int _order;
  private String _microarchitecture;
  private String _pathToLibxsmm;

  public OptimisedFluxesNonlinearADER_DGinC(int dimensions, int numberOfUnknowns, int numberOfParameters, int order,
      String microarchitecture, String pathToLibxsmm) {
    _dimensions = dimensions;
    _numberOfUnknowns = numberOfUnknowns;
    _numberOfParameters = numberOfParameters;
    _order = order;
    _microarchitecture = microarchitecture;
    _pathToLibxsmm = pathToLibxsmm;
  }

  public void writeHeader(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    Helpers.writeMinimalADERDGSolverHeader(solverName, writer, projectName);

    writer.write("  private:\n");
    if (_dimensions == 2) {
      writer.write("    static void flux(const double* const Q, double* f, double* g);\n");
    } else {
      writer.write(
          "    static void flux(const double* const Q, double* f, double* g, double* h);\n");
    }
    writer.write(
        "    static void eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda);\n");
    writer.write(
        "    static void adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q);\n");

    writer.write("};\n\n\n");
  }

  public void writeGeneratedImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException {
    Helpers.invokeCodeGenerator(solverName, _numberOfUnknowns, _numberOfParameters, _order, false, _dimensions,
        _microarchitecture, _pathToLibxsmm);

    writer.write("// ==============================================\n");
    writer.write("// Please do not change the implementations below\n");
    writer.write("// =============================---==============\n");
    writer.write("#include \"" + solverName + ".h\"\n");
    writer.write("#include \"kernels/aderdg/optimised/Kernels.h\"\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::spaceTimePredictor( double* lQi, double* lFi, double* lQhi, double* lFhi, double* lQhbnd, double* lFhbnd, const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx, const double dt ) {\n");
    writer.write("   kernels::aderdg::optimised::picardLoop<flux>( lQi, lFi, luh, dx, dt );\n");
    writer.write("   kernels::aderdg::optimised::predictor( lQhi, lFhi, lQi, lFi );\n");
    writer.write("   kernels::aderdg::optimised::extrapolator( lQhbnd, lFhbnd, lQhi, lFhi );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::solutionUpdate(double* luh, const double* const lduh, const double dt) {\n");
    writer.write("   kernels::aderdg::optimised::solutionUpdate( luh, lduh, dt );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::volumeIntegral(double* lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    writer.write("   kernels::aderdg::optimised::volumeIntegral( lduh, lFhi, dx );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::surfaceIntegral(double* lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    writer.write("   kernels::aderdg::optimised::surfaceIntegral( lduh, lFhbnd, dx );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::riemannSolver(double* FL, double* FR, const double* const QL, const double* const QR, const double dt, const int normalNonZeroIndex) {\n");
    writer.write(
        "   kernels::aderdg::optimised::riemannSolver<eigenvalues>( FL, FR, QL, QR, dt, normalNonZeroIndex );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("double " + projectName + "::" + solverName
        + "::stableTimeStepSize(const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    writer.write(
        "   return kernels::aderdg::optimised::stableTimeStepSize<eigenvalues>( luh, dx );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::solutionAdjustment(double *luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {\n");
    writer.write(
        "   kernels::aderdg::optimised::solutionAdjustment<adjustedSolutionValues>( luh, center, dx, t, dt );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
            + "::faceUnknownsProlongation(double* lQhbndFine,double* lFhbndFine,const double* lQhbndCoarse,const double* lFhbndCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) {\n");
    writer.write("   kernels::aderdg::optimised::faceUnknownsProlongation( lQhbndFine, lFhbndFine, lQhbndCoarse, lFhbndCoarse, coarseGridLevel, fineGridLevel, subfaceIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
            + "::faceUnknownsRestriction(double* lQhbndCoarse,double* lFhbndCoarse,const double* lQhbndFine,const double* lFhbndFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) {\n");
    writer.write("   kernels::aderdg::optimised::faceUnknownsRestriction( lQhbndCoarse, lFhbndCoarse, lQhbndFine, lFhbndFine, coarseGridLevel, fineGridLevel, subfaceIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
            + "::volumeUnknownsProlongation(double* luhFine, const double* luhCoarse, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {\n");
    writer.write("   kernels::aderdg::optimised::volumeUnknownsProlongation( luhFine, luhCoarse, coarseGridLevel, fineGridLevel, subcellIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
            + "::volumeUnknownsRestriction(double* luhCoarse, const double* luhFine, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {\n");
    writer.write("   kernels::aderdg::optimised::volumeUnknownsRestriction( luhCoarse, luhFine, coarseGridLevel, fineGridLevel, subcellIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
  }

  public void writeUserImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException {
    Helpers.writeMinimalADERDGSolverUserImplementation(
        solverName, writer, projectName, _numberOfUnknowns, _numberOfParameters, _order);

    int digits = String.valueOf(_numberOfUnknowns + _numberOfParameters).length();

    if (_dimensions == 2) {
      writer.write("void " + projectName + "::" + solverName
          + "::flux(const double* const Q, double* f, double* g) {\n");
    } else {
      writer.write("void " + projectName + "::" + solverName
          + "::flux(const double* const Q, double* f, double* g, double* h) {\n");
    }
    writer.write("  // Dimensions             = " + _dimensions + "\n");
    writer.write("  // Number of variables    = " + Integer.toString(_numberOfUnknowns, _numberOfParameters) + " (unknowns + parameters) \n");
    writer.write("  // f\n");
    writer.write("  // @todo Please implement\n");
    for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
      writer.write("  f[" + String.format("%" + digits + "d", i) + "] = 0.0;\n");
    }
    writer.write("  // g\n");
    writer.write("  // @todo Please implement\n");
    for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
      writer.write("  g[" + String.format("%" + digits + "d", i) + "] = 0.0;\n");
    }
    if (_dimensions == 3) {
      writer.write("  // h\n");
      writer.write("  // @todo Please implement\n");
      for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
        writer.write("  h[" + String.format("%" + digits + "d", i) + "] = 0.0;\n");
      }
    }
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {\n");
    writer.write("  // Dimensions             = " + _dimensions + "\n");
    writer.write("  // Number of variables    = " + Integer.toString(_numberOfUnknowns, _numberOfParameters) + " (unknowns + parameters) \n");
    writer.write("  // @todo Please implement\n");
    for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
      writer.write("  lambda[" + String.format("%" + digits + "d", i) + "] = 0.0;\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {\n");
    writer.write("  // Dimensions             = " + _dimensions + "\n");
    writer.write("  // Number of variables    = " + Integer.toString(_numberOfUnknowns, _numberOfParameters) + " (unknowns + parameters) \n");
    writer.write("  // @todo Please implement\n");
    for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
      writer.write("  Q[" + String.format("%" + digits + "d", i) + "] = 0.0;\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
  }
  public void writeUserPDE(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    // @todo Implement
    System.err.println("C-style kernels do not have a PDF.f90.\n");
  }
  public void writeTypesDef(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    // @todo Implement
    System.err.println("C-style kernels do not have a typesDef.f90.\n");
  }
}
