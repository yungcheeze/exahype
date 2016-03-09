package eu.exahype.solvers;

public class GenericFluxesNonlinearADER_DGinC implements Solver {
  public static final String Identifier = "generic::fluxes::nonlinear";
	  
  private int _dimensions;
  private int _numberOfVariables;
  private int _order;
  
  public GenericFluxesNonlinearADER_DGinC(int dimensions, int numberOfVariables, int order) {
	_dimensions        = dimensions;  
	_numberOfVariables = numberOfVariables;
	_order             = order;
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
	        writer.write("   kernels::aderdg::generic::c::spaceTimePredictor<flux>( lQi, lFi, lQhi, lFhi, lQhbnd, lFhbnd, luh, dx, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
	        writer.write("}\n");
	        writer.write("\n\n\n");
	        writer.write( "void " + projectName + "::" + solverName + "::solutionUpdate(double* luh, const double* const lduh, const double dt) {\n");
	        writer.write("   kernels::aderdg::generic::c::solutionUpdate( luh, lduh, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
	        writer.write("}\n");
	        writer.write("\n\n\n");
	        writer.write( "void " + projectName + "::" + solverName + "::volumeIntegral(double* lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
	        writer.write("   kernels::aderdg::generic::c::volumeIntegral( lduh, lFhi, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
	        writer.write("}\n");
	        writer.write("\n\n\n");
	        writer.write( "void " + projectName + "::" + solverName + "::surfaceIntegral(double* lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
	        writer.write("   kernels::aderdg::generic::c::surfaceIntegral( lduh, lFhbnd, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
	        writer.write("}\n");
	        writer.write("\n\n\n");
	        writer.write( "void " + projectName + "::" + solverName + "::riemannSolver(double* FL, double* FR, const double* const QL, const double* const QR, const double dt, const int normalNonZeroIndex) {\n");
	        writer.write("   kernels::aderdg::generic::c::riemannSolver<eigenvalues>( FL, FR, QL, QR, dt, normalNonZeroIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
	        writer.write("}\n");
	        writer.write("\n\n\n");
	        writer.write( "double " + projectName + "::" + solverName + "::stableTimeStepSize(const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
	        writer.write("   return kernels::aderdg::generic::c::stableTimeStepSize<eigenvalues>( luh, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
	        writer.write("}\n");
	        writer.write("\n\n\n");
	        writer.write( "void " + projectName + "::" + solverName + "::solutionAdjustment(double *luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {\n");
	        writer.write("   kernels::aderdg::generic::c::solutionAdjustment<adjustedSolutionValues>( luh, center, dx, t, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
	        writer.write("}\n");
	        writer.write("\n\n\n");
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
