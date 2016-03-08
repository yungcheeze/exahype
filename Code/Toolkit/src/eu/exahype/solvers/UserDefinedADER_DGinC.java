package eu.exahype.solvers;

public class UserDefinedADER_DGinC implements Solver {
  public static final String Identifier = "user::defined";
  
  private int _numberOfVariables;
  private int _order;
  
  public UserDefinedADER_DGinC(int numberOfVariables, int order) {
	_numberOfVariables = numberOfVariables;
	_order             = order;
  }
  
  public void writeHeader( java.io.BufferedWriter writer, String solverName, String projectName ) throws java.io.IOException {
    Helpers.writeMinimalADERDGSolverHeader( solverName, writer, projectName );

    writer.write("};\n\n\n");
  }

  public void writeGeneratedImplementation( java.io.BufferedWriter writer, String solverName, String projectName ) throws java.io.IOException {
    writer.write("// ==============================================\n");
    writer.write("// Please do not change the implementations below\n");
    writer.write("// =============================---==============\n");
    writer.write("#include \"" + solverName + ".h\"\n");
    writer.write("\n\n\n");
    writer.write("// This file is empty as a user::defined kernel is chosen, i.e. the user\n");
    writer.write("// wants to implement everything.");
    writer.write("\n\n\n");
  }
  
  public void writeUserImplementation( java.io.BufferedWriter writer, String solverName, String projectName ) throws java.io.IOException {
      writer.write("#include \"" + solverName + ".h\"\n");
      writer.write("\n\n\n");
      writer.write(projectName + "::" + solverName + "::" + solverName + "( int kernelNumber):\n");
      writer.write("  exahype::solvers::Solver(\"" + solverName + "\",exahype::solvers::Solver::ADER_DG,kernelNumber," + _numberOfVariables + "," + _order + "+1,exahype::solvers::Solver::GlobalTimeStepping) {\n");
      writer.write("  // @todo Please implement/augment if required\n");
      writer.write("}\n");
      writer.write("\n\n\n");
      writer.write("int " + projectName + "::" + solverName + "::getMinimumTreeDepth() const {\n");
      writer.write("  // @todo Please implement\n");
      writer.write("  return 3;\n");
      writer.write("}\n");
      writer.write("\n\n\n");

          writer.write("void " + projectName + "::" + solverName + "::spaceTimePredictor(double* lQi, double* lFi, double* lQhi, double* lFhi, double* lQhbnd, double* lFhbnd, const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx, const double dt ) {\n");
          writer.write("  // @todo Please implement\n");
          writer.write("}\n");
          writer.write("\n\n\n");
          writer.write("void " + projectName + "::" + solverName + "::solutionUpdate(double* luh, const double* const lduh, const double dt) {\n");
          writer.write("  // @todo Please implement\n");
          writer.write("}\n");
          writer.write("\n\n\n");
          writer.write("void " + projectName + "::" + solverName + "::volumeIntegral(double* lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
          writer.write("  // @todo Please implement\n");
          writer.write("}\n");
          writer.write("\n\n\n");
          writer.write("void " + projectName + "::" + solverName + "::surfaceIntegral(double* lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
          writer.write("  // @todo Please implement\n");
          writer.write("}\n");
          writer.write("\n\n\n");
          writer.write("void " + projectName + "::" + solverName + "::riemannSolver(double* FL, double* FR, const double* const QL, const double* const QR, const double dt, const int normalNonZeroIndex) {\n");
          writer.write("  // @todo Please implement\n");
          writer.write("}\n");
          writer.write("\n\n\n");
          writer.write("double " + projectName + "::" + solverName + "::stableTimeStepSize(const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx ) {\n");
          writer.write("  // @todo Please implement\n");
          writer.write("  return 1.0;\n");
          writer.write("}\n");
          writer.write("\n\n\n");
          writer.write("void " + projectName + "::" + solverName + "::solutionAdjustment(double* luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {\n");
          writer.write("  // @todo Please implement\n");
          writer.write("}\n");
          writer.write("\n\n\n");
  }
}
