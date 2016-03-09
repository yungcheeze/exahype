package eu.exahype.solvers;

public class GenericFluxesLinearADER_DGinC implements Solver {
  public static final String Identifier = "generic::fluxes::linear";
	  
  private int _dimensions;
  
  public GenericFluxesLinearADER_DGinC(int dimensions) {
	_dimensions = dimensions;  
  }

  public void writeHeader( java.io.BufferedWriter writer, String solverName, String projectName ) throws java.io.IOException {
    Helpers.writeMinimalADERDGSolverHeader( solverName, writer, projectName );

    writer.write("  private:\n");
    if (_dimensions==2) {
      writer.write("    static void flux(const double* const Q, double* f, double* g);\n");
    }
    else {
     writer.write("    static void flux(const double* const Q, double* f, double* g, double* h);\n");
    }
    writer.write("    static void eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda);\n");
    writer.write("    static void adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q);\n" );

    writer.write("};\n\n\n");
  }

	  public void writeGeneratedImplementation( java.io.BufferedWriter writer, String solverName, String projectName ) throws java.io.IOException {
	    // @todo Implement
		System.err.println( "not implemented yet\n" );
	  }
	  
	  public void writeUserImplementation( java.io.BufferedWriter writer, String solverName, String projectName ) throws java.io.IOException {
	    // @todo Implement
	    System.err.println( "not implemented yet\n" );
	  }
}
