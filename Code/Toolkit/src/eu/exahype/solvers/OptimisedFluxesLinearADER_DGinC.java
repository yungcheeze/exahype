package eu.exahype.solvers;

public class OptimisedFluxesLinearADER_DGinC implements Solver {
  public static final String Identifier = "optimised::fluxes::linear";
	  
  private int _dimensions;
  
  public OptimisedFluxesLinearADER_DGinC(int dimensions) {
	_dimensions = dimensions;  
  }

  public void writeHeader( java.io.BufferedWriter writer, String solverName, String projectName ) throws java.io.IOException {
		// @todo
		System.err.println( "not implemented yet\n" );
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
