package eu.exahype.solvers;

public class GenericFluxesNonlinearADER_DGinFortran implements Solver {
	  public static final String Identifier = GenericFluxesNonlinearADER_DGinC.Identifier;
	  
  private int _dimensions;
  private int _numberOfVariables;
  private int _order;
  
  public GenericFluxesNonlinearADER_DGinFortran(int dimensions, int numberOfVariables, int order) {
	_dimensions        = dimensions;  
	_numberOfVariables = numberOfVariables;
	_order             = order;
  }
	  public void writeHeader( java.io.BufferedWriter writer, String solverName, String projectName ) throws java.io.IOException {
		// @todo
    Helpers.writeMinimalADERDGSolverHeader( solverName, writer, projectName );

    writer.write("  private:\n");
    if (_dimensions==2) {
      writer.write("    static void flux(const double* const Q, double* f, double* g);\n");
    }
    else {
      writer.write("    static void flux(const double* const Q, double* f, double* g, double* h);\n");
    }
    writer.write("    static void eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda);\n");
    writer.write("    static void adjustedSolutionValues(const double* const x,const double J_w,const double t,const double dt,double* Q);\n" );

    writer.write("};\n\n\n");
  }

	  public void writeGeneratedImplementation( java.io.BufferedWriter writer, String solverName, String projectName ) throws java.io.IOException {
	        writer.write("// ==============================================\n");
	        writer.write("// Please do not change the implementations below\n");
	        writer.write("// =============================---==============\n");
	        writer.write("#include \"" + solverName + ".h\"\n");
	        writer.write("#include \"kernels/aderdg/generic/Kernels.h\"\n");
	        writer.write("\n\n\n");
	        writer.write( "void " + projectName + "::" + solverName + "::spaceTimePredictor( double* lQi, double* lFi, double* lQhi, double* lFhi, double* lQhbnd, double* lFhbnd, const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx, const double dt ) {\n");
	        writer.write("   kernels::aderdg::generic::fortran::spaceTimePredictor<flux>( lQi, lFi, lQhi, lFhi, lQhbnd, lFhbnd, luh, dx, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
	        writer.write("}\n");
	        writer.write("\n\n\n");
	        writer.write( "void " + projectName + "::" + solverName + "::solutionUpdate(double* luh, const double* const lduh, const double dt) {\n");
	        writer.write("   kernels::aderdg::generic::fortran::solutionUpdate( luh, lduh, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
	        writer.write("}\n");
	        writer.write("\n\n\n");
	        writer.write( "void " + projectName + "::" + solverName + "::volumeIntegral(double* lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
	        writer.write("   kernels::aderdg::generic::fortran::volumeIntegral( lduh, lFhi, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
	        writer.write("}\n");
	        writer.write("\n\n\n");
	        writer.write( "void " + projectName + "::" + solverName + "::surfaceIntegral(double* lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
	        writer.write("   kernels::aderdg::generic::fortran::surfaceIntegral( lduh, lFhbnd, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
	        writer.write("}\n");
	        writer.write("\n\n\n");
	        writer.write( "void " + projectName + "::" + solverName + "::riemannSolver(double* FL, double* FR, const double* const QL, const double* const QR, const double dt, const int normalNonZeroIndex) {\n");
	        writer.write("   kernels::aderdg::generic::fortran::riemannSolver<eigenvalues>( FL, FR, QL, QR, dt, normalNonZeroIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
	        writer.write("}\n");
	        writer.write("\n\n\n");
	        writer.write( "double " + projectName + "::" + solverName + "::stableTimeStepSize(const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
	        writer.write("   return kernels::aderdg::generic::fortran::stableTimeStepSize<eigenvalues>( luh, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
	        writer.write("}\n");
	        writer.write("\n\n\n");
	        writer.write( "void " + projectName + "::" + solverName + "::solutionAdjustment(double *luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {\n");
	        writer.write("   kernels::aderdg::generic::fortran::solutionAdjustment<adjustedSolutionValues>( luh, center, dx, t, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
	        writer.write("}\n");
	        writer.write("\n\n\n");
	  }
	  
	  public void writeUserImplementation( java.io.BufferedWriter writer, String solverName, String projectName ) throws java.io.IOException {
	    // @todo Implement
	    System.err.println( "not implemented yet\n" );
	  }
}
