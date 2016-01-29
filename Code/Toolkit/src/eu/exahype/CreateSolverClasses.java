package eu.exahype;


import java.io.IOException;

import eu.exahype.analysis.DepthFirstAdapter;
import eu.exahype.node.AAderdgSolver;
import eu.exahype.node.AProject;
import eu.exahype.node.ATwoDimensionalComputationalDomain;
import eu.exahype.node.AThreeDimensionalComputationalDomain;


public class CreateSolverClasses extends DepthFirstAdapter {
  public Boolean valid = true;

  private DirectoryAndPathChecker   _directoryAndPathChecker;
  
  private String                    _projectName;
  
  private int                       _dimensions;
  

  public CreateSolverClasses(DirectoryAndPathChecker  directoryAndPathChecker) {
	_directoryAndPathChecker = directoryAndPathChecker;
  }

  @Override
  public void inAProject(AProject node) {
	_projectName = node.getName().toString().trim();
	
    if (node.getSolver().size()==0) { 
      System.out.println( "there are no solvers in the specification file ... nothing to be done" );      
    }
  }
  

  @Override
  public void inATwoDimensionalComputationalDomain(ATwoDimensionalComputationalDomain node) {
    _dimensions = 2;
  }
  
  
  @Override
  public void inAThreeDimensionalComputationalDomain(AThreeDimensionalComputationalDomain node) {
    _dimensions = 3;
  }
  
    
  @Override
  public void inAAderdgSolver(AAderdgSolver node) {
	String solverName = node.getName().toString().trim();

	java.io.File headerFile                  = new java.io.File(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/" + solverName + ".h");
	java.io.File userImplementationFile      = new java.io.File(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/" + solverName + ".cpp");
	java.io.File generatedImplementationFile = new java.io.File(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/" + solverName + "_generated.cpp");
	
	String kernel = node.getKernel().toString().trim();
	
	try {	  
      if (headerFile.exists()) {
        System.out.println( "create header of solver " + solverName + " ... header " + headerFile.getAbsoluteFile() + " does exist already. Remove to allow toolkit to regenerate it (changes will be lost)" );      
	  }
	  else {
		java.io.BufferedWriter headerWriter = new java.io.BufferedWriter(new java.io.FileWriter(headerFile));
		
		if (kernel.equals("user::defined")) {
          System.out.println( "create header of solver " + solverName + " ... user::defined kernel" );      

          writeMinimalADERDGSolverHeader( solverName, headerWriter );
        }
		else if (kernel.equals("user::fluxes")) {
          System.out.println( "create header of solver " + solverName + " ... user::fluxes kernel" );      

          writeADERDGSolverHeaderForUserFluxes( solverName, headerWriter );
		}
		// Vasco, Angelika: your turn
/*		else if (node.getKernel().toString().trim().equals("user::fluxes")) {
          System.out.println( "create header of solver " + solverName + " ... user::fluxes kernel" );      
          writeMinimalADERDGSolverHeader( solverName, headerWriter );
		}*/
		else {
	      System.err.println( "ERROR: unknown ADER-DG kernel type " + kernel + " ... user::defined kernel" );      
		  valid = false;
		}
		
		headerWriter.close();
	  }
	
      if (userImplementationFile.exists()) {
        System.out.println( "user's implementation file of solver " + solverName + " ... does exist already. Is not overwritten" );      
      }
	  else {
        java.io.BufferedWriter userImplementationWriter = new java.io.BufferedWriter(new java.io.FileWriter(userImplementationFile));
		
 		if (kernel.equals("user::defined")) {
           System.out.println( "create header of solver " + solverName + " ... user::defined kernel" );      
           
           writeADERDGSolverUserImplementationForUserDefined(solverName, node.getVariables().toString().trim(), node.getOrder().toString().trim(), userImplementationWriter);
 		}
 		else if (kernel.equals("user::fluxes")) {
           System.out.println( "create header of solver " + solverName + " ... user::fluxes kernel" );      

           writeADERDGSolverUserImplementationForUserFluxes(solverName, node.getVariables().toString().trim(), node.getOrder().toString().trim(), userImplementationWriter);
 		}
 		// Vasco, Angelika: your turn
 /*		else if (node.getKernel().toString().trim().equals("user::fluxes")) {
           System.out.println( "create header of solver " + solverName + " ... user::fluxes kernel" );      
 		}*/
 		else {
 	      System.err.println( "ERROR: unknown ADER-DG kernel type " + kernel + " ... user::defined kernel" );      
 		  valid = false;
 		}

 		userImplementationWriter.close();
	  }

	  if (generatedImplementationFile.exists()) {
	    System.out.println( "generated implementation file of solver " + solverName + " ... does exist already. Is overwritten" );      
	  }
	  java.io.BufferedWriter generatedImplementationWriter = new java.io.BufferedWriter(new java.io.FileWriter(generatedImplementationFile));

      if (kernel.equals("user::defined")) {
        System.out.println( "create header of solver " + solverName + " ... user::defined kernel" );      
        
        writeADERDGSolverGeneratedImplementationForUserDefined(solverName, generatedImplementationWriter);
	  }
      else if (kernel.equals("user::fluxes")) {
        System.out.println( "create header of solver " + solverName + " ... user::fluxes kernel" );      

        writeADERDGSolverGeneratedImplementationForUserFluxes(solverName, node.getVariables().toString().trim(), node.getOrder().toString().trim(), generatedImplementationWriter);
      }
		// Vasco, Angelika: your turn
/*		else if (node.getKernel().toString().trim().equals("user::fluxes")) {
        System.out.println( "create header of solver " + solverName + " ... user::fluxes kernel" );      
		}*/
	  else {
	    System.err.println( "ERROR: unknown ADER-DG kernel type " + kernel + " ... user::defined kernel" );      
        valid = false;
	  }
      
      generatedImplementationWriter.close();
	}	  
	catch (Exception exc) {
      System.err.println( "ERROR: " + exc.toString() );
	  valid = false;
	}
  }


  private void writeMinimalADERDGSolverHeader(
	String                 solverName,
	java.io.BufferedWriter writer
  ) throws IOException {
	  writer.write("// This file is generated by the ExaHyPE toolkit.\n");
	  writer.write("// Please do not modify - it will be overwritten by the next\n");
	  writer.write("// ExaHyPE toolkit call.\n");
	  writer.write("// \n");
	  writer.write("// ========================\n");
	  writer.write("//   www.exahype.eu\n");
	  writer.write("// ========================\n");

	  writer.write("\n\n\n");
	  writer.write("#include \"exahype/solvers/Solver.h\"");
      writer.write("\n\n\n");
	  
	  writer.write("namespace " + _projectName + "{\n");
	  writer.write("  class " + solverName + ";\n");
	  writer.write("}\n\n\n");
	  
	  writer.write("class " + _projectName + "::" + solverName + ": public exahype::solvers::Solver {\n");
	  writer.write("  public:\n");
	  writer.write("    " + solverName + "(int kernelNumber); \n");
	  writer.write("    virtual int getMinimumTreeDepth() const;\n");
	  
	  if (_dimensions==2) {
        writer.write("    virtual void spaceTimePredictor( double * lQi, double * lFi, const double * const luh, double * lQhi, double * lFhi, double * lQhbnd, double * lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& x, const double dt ); \n");
	  }
	  else {
        writer.write("ERROR\n");
	  }
	  
	  writer.write("};\n\n\n");
  }


  private void writeADERDGSolverHeaderForUserFluxes(
	String                 solverName,
	java.io.BufferedWriter writer
  ) throws IOException {
	  writer.write("// This file is generated by the ExaHyPE toolkit.\n");
	  writer.write("// Please do not modify - it will be overwritten by the next\n");
	  writer.write("// ExaHyPE toolkit call.\n");
	  writer.write("// \n");
	  writer.write("// ========================\n");
	  writer.write("//   www.exahype.eu\n");
	  writer.write("// ========================\n");

	  writer.write("\n\n\n");
	  writer.write("#include \"exahype/solvers/Solver.h\"");
      writer.write("\n\n\n");
	  
	  writer.write("namespace " + _projectName + "{\n");
	  writer.write("  class " + solverName + ";\n");
	  writer.write("}\n\n\n");
	  
	  writer.write("class " + _projectName + "::" + solverName + ": public exahype::solvers::Solver {\n");
	  writer.write("  private:\n");
      if (_dimensions==2) {
	    writer.write("    static void flux(const double * const Q, double * f, double * g);\n");
	    writer.write("    static void eigenvalues(const double * const Q, const int normalNonZeroIndex, double * lambda);\n");
	  }
	  else {
	    writer.write("ERROR\n");
	  }
	  writer.write("  public:\n");
	  writer.write("    " + solverName + "(int kernelNumber); \n");
	  writer.write("    virtual int getMinimumTreeDepth() const;\n");
	  
	  if (_dimensions==2) {
        writer.write("    virtual void spaceTimePredictor(double * lQi, double * lFi, double * lQhi, double * lFhi, double * lQhbnd, double * lFhbnd, const double * const luh, const tarch::la::Vector<DIMENSIONS,double>& dx, const double dt ); \n");
        writer.write("    virtual void solutionUpdate(double * luh, const double * const lduh, const tarch::la::Vector<DIMENSIONS,double>& dx, const double dt);\n");
        writer.write("    virtual void volumeIntegral(double * lduh, const double * const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx);\n" );
        writer.write("    virtual void surfaceIntegral(double * lduh, const double * const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx);\n" );
        writer.write("    virtual void riemannSolver(double * FL, double * FR, const double * const QL, const double * const QR, const double dt, const int normalNonZeroIndex);\n" );
        writer.write("    virtual double stableTimeStepSize(const double * const luh, const tarch::la::Vector<DIMENSIONS,double>& dx );\n" );
        writer.write("    virtual void initialValues(double * luh, const tarch::la::Vector<DIMENSIONS,double>&  center, const tarch::la::Vector<DIMENSIONS,double>& dx);\n" );
	  }
	  else {
        writer.write("ERROR\n");
	  }
	  
	  writer.write("};\n\n\n");
  }
  
    
  private void writeADERDGSolverGeneratedImplementationForUserDefined(
	String                 solverName,
	java.io.BufferedWriter writer
  ) throws IOException {
    writer.write("// ==============================================\n");
  	writer.write("// Please do not change the implementations below\n");
	writer.write("// =============================---==============\n");
	writer.write("#include \"" + solverName + ".h\"\n");
    writer.write("\n\n\n");
  }

  
  private void writeADERDGSolverGeneratedImplementationForUserFluxes(
	String                 solverName,
	String                 numberOfVariables,
	String                 order,
	java.io.BufferedWriter writer
  ) throws IOException {
    writer.write("// ==============================================\n");
  	writer.write("// Please do not change the implementations below\n");
	writer.write("// =============================---==============\n");
	writer.write("#include \"" + solverName + ".h\"\n");
	writer.write("#include \"kernels/aderdg/generic/Kernels.h\"\n");
    writer.write("\n\n\n");
    writer.write( "#include \"kernels/aderdg/generic/Kernels.h\"\n\n\n");
    writer.write( "void " + _projectName + "::" + solverName + "::spaceTimePredictor( double * lQi, double * lFi, double * lQhi, double * lFhi, double * lQhbnd, double * lFhbnd, const double * const luh, const tarch::la::Vector<DIMENSIONS,double>& dx, const double dt ) {\n");
    writer.write("   kernels::aderdg::generic::spaceTimePredictor<flux>( lQi,lFi,luh,lQhi,lFhi,lQhbnd,lFhbnd,dx,dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write( "void " + _projectName + "::" + solverName + "::solutionUpdate(double * luh, const double * const lduh, const tarch::la::Vector<DIMENSIONS,double>& dx, const double dt) {\n");
    writer.write("   kernels::aderdg::generic::solutionUpdate( luh, lduh,dx, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write( "void " + _projectName + "::" + solverName + "::volumeIntegral(double * lduh, const double * const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    writer.write("   kernels::aderdg::generic::volumeIntegral( lduh, lFhi,dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write( "void " + _projectName + "::" + solverName + "::surfaceIntegral(double * lduh, const double * const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    writer.write("   kernels::aderdg::generic::surfaceIntegral( lduh, lFhbnd, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write( "void " + _projectName + "::" + solverName + "::riemannSolver(double * FL, double * FR, const double * const QL, const double * const QR, const double dt, const int normalNonZeroIndex) {\n");
    writer.write("   kernels::aderdg::generic::riemannSolver<eigenvalues>( FL, FR, QL, QR, dt, normalNonZeroIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write( "double " + _projectName + "::" + solverName + "::stableTimeStepSize(const double * const luh, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    writer.write("   return kernels::aderdg::generic::stableTimeStepSize<eigenvalues>( luh, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    writer.write("}\n");
    writer.write("\n\n\n");
  }

  
  private void writeADERDGSolverUserImplementationForUserDefined(
	String                 solverName,
	String                 numberOfVariables,
	String                 order,
	java.io.BufferedWriter writer
  ) throws IOException {
	writer.write("#include \"" + solverName + ".h\"\n");
    writer.write("\n\n\n");
    writer.write(_projectName + "::" + solverName + "::" + solverName + "( int kernelNumber):\n");
    writer.write("  exahype::solvers::Solver(\"" + solverName + "\",exahype::solvers::Solver::ADER_DG,kernelNumber," + numberOfVariables + "," + order + "+1) {\n");
    writer.write("  // @todo Please implement/augment if required\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("int " + _projectName + "::" + solverName + "::getMinimumTreeDepth() const {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("  return 3;\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write( "#include \"kernels/aderdg/generic/Kernels.h\"\n\n\n");
    writer.write( "void " + _projectName + "::" + solverName + "::spaceTimePredictor( double * lQi, double * lFi, const double * const luh, double * lQhi, double * lFhi, double * lQhbnd, double * lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& x, const double dt ) {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("}\n");
    writer.write("\n\n\n");
  }

  
  private void writeADERDGSolverUserImplementationForUserFluxes(
	String                 solverName,
	String                 numberOfVariables,
	String                 order,
	java.io.BufferedWriter writer
  ) throws IOException {
	writer.write("#include \"" + solverName + ".h\"\n");
    writer.write("\n\n\n");
    writer.write(_projectName + "::" + solverName + "::" + solverName + "( int kernelNumber):\n");
    writer.write("  exahype::solvers::Solver(\"" + solverName + "\",exahype::solvers::Solver::ADER_DG,kernelNumber," + numberOfVariables + "," + order + "+1) {\n");
    writer.write("  // @todo Please implement/augment if required\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("int " + _projectName + "::" + solverName + "::getMinimumTreeDepth() const {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("  return 3;\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + _projectName + "::" + solverName + "::flux(const double * const Q, double * f, double * g) {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + _projectName + "::" + solverName + "::eigenvalues(const double * const Q, const int normalNonZeroIndex, double * lambda) {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + _projectName + "::" + solverName + "::initialValues(double * luh, const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("}\n");
    writer.write("\n\n\n");
  }
} 
