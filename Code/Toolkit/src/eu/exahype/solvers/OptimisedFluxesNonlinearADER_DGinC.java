package eu.exahype.solvers;

public class OptimisedFluxesNonlinearADER_DGinC implements Solver {
  public static final String Identifier = "optimised::fluxes::nonlinear";

  private int _dimensions;
  private int _numberOfVariables;
  private int _order;
  
  public OptimisedFluxesNonlinearADER_DGinC(int dimensions, int numberOfVariables, int order) {
	_dimensions        = dimensions;  
	_numberOfVariables = numberOfVariables;
	_order             = order;
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
	        Helpers.writeMinimalADERDGSolverUserImplementation(solverName,writer, projectName, _numberOfVariables,_order );

          int digits = String.valueOf(_numberOfVariables).length();

          if (_dimensions==2) {
            writer.write("void " + projectName + "::" + solverName + "::flux(const double* const Q, double* f, double* g) {\n");
          }
          else {
            writer.write("void " + projectName + "::" + solverName + "::flux(const double* const Q, double* f, double* g, double* h) {\n");
          }
          writer.write("  // Dimensions             = "+_dimensions      +"\n");
          writer.write("  // Number of variables    = "+Integer.toString(_numberOfVariables)+"\n");
          writer.write("  // f\n");
          writer.write("  // @todo Please implement\n");
          for (int i=0; i < _numberOfVariables; i++) {
              writer.write("  f["+String.format("%"+digits+"d",i)+"] = 0.0;\n");    
          }
          writer.write("  // g\n");
          writer.write("  // @todo Please implement\n");
          for (int i=0; i < _numberOfVariables; i++) {
              writer.write("  g["+String.format("%"+digits+"d",i)+"] = 0.0;\n");    
          }
          if (_dimensions==3) {
            writer.write("  // h\n");
            writer.write("  // @todo Please implement\n");
            for (int i=0; i < _numberOfVariables; i++) {
                writer.write("  h["+String.format("%"+digits+"d",i)+"] = 0.0;\n");    
            }
          }
          writer.write("}\n");
          writer.write("\n\n\n");
          writer.write("void " + projectName + "::" + solverName + "::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {\n");
          writer.write("  // Dimensions             = "+_dimensions      +"\n");
          writer.write("  // Number of variables    = "+Integer.toString(_numberOfVariables)+"\n");
          writer.write("  // @todo Please implement\n");
          for (int i=0; i < _numberOfVariables; i++) {
              writer.write("  lambda["+String.format("%"+digits+"d",i)+"] = 0.0;\n");    
          }
          writer.write("}\n");
          writer.write("\n\n\n");
          writer.write("void " + projectName + "::" + solverName + "::adjustedSolutionValues(const double* const x,const double J_w,const double t,const double dt,double* Q) {\n");
          writer.write("  // Dimensions             = "+_dimensions      +"\n");
          writer.write("  // Number of variables    = "+Integer.toString(_numberOfVariables)+"\n");
          writer.write("  // @todo Please implement\n");
          for (int i=0; i < _numberOfVariables; i++) {
              writer.write("  Q["+String.format("%"+digits+"d",i)+"] = 0.0;\n");    
          }
          writer.write("}\n");
          writer.write("\n\n\n");
	  }
}
