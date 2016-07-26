package eu.exahype.solvers;

public class GenericFiniteVolumesMUSCLinC implements Solver {
  public static final String Identifier = "generic::MUSCL";

  private int _numberOfVariables;
  private int _numberOfParameters;
  private int _patchSize;

  public GenericFiniteVolumesMUSCLinC(int numberOfVariables, int numberOfParameters, int patchSize) {
    _numberOfVariables  = numberOfVariables;
    _numberOfParameters = numberOfParameters;
    _patchSize = patchSize;
  }

  public void writeHeader(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    Helpers.writeMinimalFiniteVolumesSolverHeader(solverName, writer, projectName);

    writer.write("  private:\n");
    writer.write("    static void eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda);\n");
    writer.write("    static void flux(const double* const Q, double** F);\n");

    writer.write("};\n\n\n");
  }

  public void writeGeneratedImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException {
    writer.write("// ==============================================\n");
    writer.write("// Please do not change the implementations below\n");
    writer.write("// =============================---==============\n");
    writer.write("#include \"" + solverName + ".h\"\n");
    writer.write("\n\n\n");

    writer.write("double " + projectName + "::" + solverName + "::stableTimeStepSize( const double* const luh, const tarch::la::Vector<DIMENSIONS, double>& dx) {\n" );
    writer.write("  assertionMsg( false, \"not yet inserted in GenericFiniteVolumesMUSCLinC\");\n" );
    writer.write("}\n\n\n" );
    
    writer.write("void " + projectName + "::" + solverName + "::solutionUpdate(double** luh, const tarch::la::Vector<DIMENSIONS, double>& dx, const double dt, double& maxAdmissibleDt) {\n" );
    writer.write("  assertionMsg( false, \"not yet inserted in GenericFiniteVolumesMUSCLinC\");\n" );
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
    writer.write("  // @todo Please implement\n");
    writer.write("}\n\n\n" );

    writer.write("bool " + projectName + "::" + solverName + "::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t) {\n" );
    writer.write("  // @todo Please implement\n");
    writer.write("  return false; \n");
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
