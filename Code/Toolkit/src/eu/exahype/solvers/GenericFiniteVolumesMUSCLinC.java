package eu.exahype.solvers;

public class GenericFiniteVolumesMUSCLinC implements Solver {
  public static final String Identifier = "generic::MUSCL";

  private int _numberOfVariables;
  private int _numberOfParameters;
  private int _patchSize;
  private boolean _enableProfiler;

  public GenericFiniteVolumesMUSCLinC(int numberOfVariables, int numberOfParameters, int patchSize, boolean enableProfiler) {
    _numberOfVariables  = numberOfVariables;
    _numberOfParameters = numberOfParameters;
    _patchSize          = patchSize;
    _enableProfiler     = enableProfiler;
  }

  public void writeHeader(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    Helpers.writeMinimalFiniteVolumesSolverHeader(solverName, writer, projectName);

    writer.write("  private:\n");
    writer.write("    static void eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda);\n");
    writer.write("    static void flux(const double* const Q, double** F);\n");
    writer.write("    static void source(const double* const Q, double* S);\n");
    writer.write("};\n\n\n");
  }

  public void writeGeneratedImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException {
    writer.write("// ==============================================\n");
    writer.write("// Please do not change the implementations below\n");
    writer.write("// =============================---==============\n");
    writer.write("#include \"" + solverName + ".h\"\n");
    writer.write("#include \"kernels/finitevolumes/muscl/c/2d/solutionUpdate.cpph\"\n");
    writer.write("#include \"kernels/finitevolumes/muscl/c/3d/solutionUpdate.cpph\"\n");
    writer.write("\n\n\n");

    writer.write("double " + projectName + "::" + solverName + "::stableTimeStepSize( double* luh[THREE_POWER_D], const tarch::la::Vector<DIMENSIONS, double>& dx) {\n" );
    if (_enableProfiler) {
        writer.write("  _profiler->start(\"solutionUpdate\");\n");
      }
      writer.write( 
        "  double maxAdmissibleDt=0;\n"
      );
      writer.write(
        "  kernels::finitevolumes::muscl::c::solutionUpdate<flux,eigenvalues>( luh, dx, std::numeric_limits<double>::max(), getNumberOfVariables(), getNodesPerCoordinateAxis(), maxAdmissibleDt );\n"
      );
      if (_enableProfiler) {
        writer.write("  _profiler->stop(\"solutionUpdate\");\n");
      }
      writer.write( 
        "  return maxAdmissibleDt;\n"
      );
    writer.write("}\n\n\n" );

    writer.write("void " + projectName + "::" + solverName + "::solutionUpdate(double* luh[THREE_POWER_D], const tarch::la::Vector<DIMENSIONS, double>& dx, const double dt, double& maxAdmissibleDt) {\n" );
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"solutionUpdate\");\n");
    }
    writer.write(
      "  kernels::finitevolumes::muscl::c::solutionUpdate<flux,eigenvalues>( luh, dx, dt, getNumberOfVariables(), getNodesPerCoordinateAxis(), maxAdmissibleDt );\n"
    );
    if (_enableProfiler) {
      writer.write("  _profiler->stop(\"solutionUpdate\");\n");
    }
    writer.write("}\n\n\n" );
    
    writer.write("\n\n\n");
  }

  public void writeUserImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException {
    writer.write("#include \"" + solverName + ".h\"\n");
    writer.write("\n\n\n");
    writer.write(projectName + "::" + solverName + "::" + solverName + "(int cellsPerCoordinateAxis, double maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping, std::unique_ptr<exahype::profilers::Profiler> profiler):\n");
    writer.write("  exahype::solvers::FiniteVolumesSolver("
            + "\""+solverName+"\", "+_numberOfVariables+", "+_numberOfParameters+", cellsPerCoordinateAxis, maximumMeshSize, timeStepping, std::move(profiler)) {\n");
    writer.write("  // @todo Please implement/augment if required\n");
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName + "::solutionAdjustment( double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, double dt) {\n" );
    writer.write("  if ( tarch::la::equals(t,0.0) ) {\n");
    writer.write("  // @todo Please implement and set initial conditions\n");
    writer.write("  }\n");
    writer.write("  // @todo Feel free to add further conditions\n");
    writer.write("}\n\n\n" );

    writer.write("bool " + projectName + "::" + solverName + "::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t) {\n" );
    writer.write("  // @todo Please implement\n");
    writer.write("  if ( tarch::la::equals(t,0.0) ) {\n");
    writer.write("    // Tell kernel that you want to set initial conditions \n");
    writer.write("    return true;\n");
    writer.write("  }\n");
    writer.write("  else {\n");
    writer.write("    // @todo Please implement\n");
    writer.write("    return false; \n");
    writer.write("  }\n");
    writer.write("}\n\n\n" );

    writer.write("exahype::solvers::Solver::RefinementControl " + projectName + "::" + solverName + "::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,const tarch::la::Vector<DIMENSIONS, double>& dx, double t,const int level) {\n" );
    writer.write("  // @todo Please implement\n");
    writer.write("  return exahype::solvers::Solver::RefinementControl::Keep;\n" );
    writer.write("}\n\n\n" );

    writer.write("void " + projectName + "::" + solverName + "::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {\n" );
    writer.write("  // @todo Please implement\n");
    writer.write("}\n\n\n" );

    writer.write("void " + projectName + "::" + solverName + "::flux(const double* const Q, double** F) {\n" );
    writer.write("  // @todo Please implement\n");
    writer.write("}\n\n\n" );
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
