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
    
    private String                    _microarchitecture;

    private java.util.List<String>    _supportedMicroarchitectures;
    
    private String                    _pathToLibxsmm;

    private int                       _dimensions;
    
    private enum                      Numerics {LINEAR, NONLINEAR} 


    public CreateSolverClasses(DirectoryAndPathChecker  directoryAndPathChecker) {
        _directoryAndPathChecker = directoryAndPathChecker;
        _supportedMicroarchitectures = java.util.Arrays.asList("wsm", "snb", "hsw", "knc", "knl", "noarch");
    }

    @Override
    public void inAProject(AProject node) {
        _projectName = node.getName().toString().trim();
        
        if (node.getSolver().size()==0) { 
            System.out.println( "there are no solvers in the specification file ... nothing to be done" );      
        }
        
        _microarchitecture = node.getArchitecture().toString().trim().toLowerCase();
        if(!_supportedMicroarchitectures.contains(_microarchitecture)) {
            System.out.println( "Unknown architecture specified ... fallback solution \"noarch\" taken" );
            _microarchitecture = "noarch";
        }

    } 
    
    
    @Override
    public void inAPaths(eu.exahype.node.APaths node) {
        if(node.getLibxsmmPath() == null) {
            // attribute 'libxsmm-path' did not occur in spec file
            _pathToLibxsmm = "";    
        } 
        else {
            _pathToLibxsmm = node.getLibxsmmPath().toString().trim();
        }
    };


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

      boolean isFortran = false;
      if ( node.getLanguage().getText().trim().equals("C") ) {
        isFortran = false;
      }
      else if ( node.getLanguage().getText().trim().equals("Fortran") ) {
       isFortran = true;		
      }
      else {
        System.err.println( "ERROR: unknown language for solver " + node.getName().getText() + ". Supported languages are C and Fortran" );
        valid = false;
        return;
      }

      int numberOfVariables = Integer.parseInt(node.getVariables().toString().trim());
      int order             = Integer.parseInt(node.getOrder().toString().trim());
      
      eu.exahype.solvers.Solver solver = null;
      if (isFortran) {
        switch(kernel) {
          case eu.exahype.solvers.UserDefinedADER_DGinFortran.Identifier:
            solver = new eu.exahype.solvers.UserDefinedADER_DGinFortran();
            break;
          case eu.exahype.solvers.GenericFluxesLinearADER_DGinFortran.Identifier:
            solver = new eu.exahype.solvers.GenericFluxesLinearADER_DGinFortran();
            break;
          case eu.exahype.solvers.GenericFluxesNonlinearADER_DGinFortran.Identifier:
            solver = new eu.exahype.solvers.GenericFluxesNonlinearADER_DGinFortran();
            break;
        }
      }
      else {
        switch(kernel) {
          case eu.exahype.solvers.UserDefinedADER_DGinC.Identifier:
            solver = new eu.exahype.solvers.UserDefinedADER_DGinC(numberOfVariables,order);
            break;
          case eu.exahype.solvers.GenericFluxesLinearADER_DGinC.Identifier:
            solver = new eu.exahype.solvers.GenericFluxesLinearADER_DGinC(_dimensions);
            break;
          case eu.exahype.solvers.GenericFluxesNonlinearADER_DGinC.Identifier:
            solver = new eu.exahype.solvers.GenericFluxesNonlinearADER_DGinC(_dimensions,numberOfVariables,order);
            break;
          case eu.exahype.solvers.OptimisedFluxesLinearADER_DGinC.Identifier:
            solver = new eu.exahype.solvers.OptimisedFluxesLinearADER_DGinC(_dimensions);
            break;
          case eu.exahype.solvers.OptimisedFluxesNonlinearADER_DGinC.Identifier:
            solver = new eu.exahype.solvers.OptimisedFluxesNonlinearADER_DGinC(_dimensions,numberOfVariables,order);
            break;
          case eu.exahype.solvers.KernelEuler2d.Identifier:
            solver = new eu.exahype.solvers.KernelEuler2d();
            break;
        }
      }
      
      if (solver==null) {
        System.err.println( "creation solver " + solverName + " ... failed as kernel " + kernel + " for language " + node.getLanguage().getText().trim() + " is not supported" );
        valid = false;
        return;
      }
      
      try {      
        // =====================
        // Write all the headers
        // =====================
        if (headerFile.exists()) {
            System.out.println( "create header of solver " + solverName + " ... header " + headerFile.getAbsoluteFile() + " does exist already. Remove to allow toolkit to regenerate it (changes will be lost)" );      
        }
        else {
          java.io.BufferedWriter headerWriter = new java.io.BufferedWriter(new java.io.FileWriter(headerFile));
          solver.writeHeader(headerWriter, solverName, _projectName);
          System.out.println( "create header of solver " + solverName + " ... ok" );
          headerWriter.close();
        }

        if (userImplementationFile.exists()) {
          System.out.println( "user's implementation file of solver " + solverName + " ... does exist already. Is not overwritten" );      
        }
        else {
          java.io.BufferedWriter userImplementationWriter = new java.io.BufferedWriter(new java.io.FileWriter(userImplementationFile));
          solver.writeUserImplementation(userImplementationWriter, solverName, _projectName);
          System.out.println( "create user implementation template of solver " + solverName + " ... please complete" );
          userImplementationWriter.close();
        }
        
        if (generatedImplementationFile.exists()) {
          System.out.println( "generated implementation file of solver " + solverName + " ... does exist already. Is overwritten" );      
        }

        java.io.BufferedWriter generatedImplementationWriter = new java.io.BufferedWriter(new java.io.FileWriter(generatedImplementationFile));
        solver.writeGeneratedImplementation(generatedImplementationWriter, solverName, _projectName);
        System.out.println( "create generated implementation of solver " + solverName + " ... ok" );
        generatedImplementationWriter.close();


/*

                case "generic::fluxes::nonlinear":
                    System.out.println( "create generated implementation of solver " + solverName + " ... ok" );
                    writeADERDGSolverGeneratedImplementationForUserFluxes(solverName, node.getVariables().toString().trim(), node.getOrder().toString().trim(), generatedImplementationWriter);
                    break;
                    
                case "optimised::fluxes::nonlinear":
                    System.out.println( "create generated implementation of solver " + solverName + " ... ok" );
                    writeADERDGSolverGeneratedImplementationForOptimisedKernel(solverName, node.getVariables().toString().trim(), node.getOrder().toString().trim(), Numerics.NONLINEAR, generatedImplementationWriter);
                    invokeCodeGenerator(solverName, node.getVariables().toString().trim(), node.getOrder().toString().trim(), Numerics.NONLINEAR);
                    break;
                    
                case "optimised::fluxes::linear":
                    System.out.println( "create generated implementation of solver " + solverName + " ... ok" );
                    writeADERDGSolverGeneratedImplementationForOptimisedKernel(solverName, node.getVariables().toString().trim(), node.getOrder().toString().trim(), Numerics.LINEAR, generatedImplementationWriter);
                    invokeCodeGenerator(solverName, node.getVariables().toString().trim(), node.getOrder().toString().trim(), Numerics.LINEAR);
                    break;
                    
                default:
                    System.err.println( "ERROR: unknown ADER-DG kernel type " + kernel + " ... failed" );      
                    valid = false;
            }

            generatedImplementationWriter.close();*/
      }      
      catch (Exception exc) {
        System.err.println( "ERROR: " + exc.toString() );
        valid = false;
      }
    }

/*

    private void writeADERDGSolverGeneratedImplementationForOptimisedKernel(
            String                 solverName,
            String                 numberOfVariables,
            String                 order,
            Numerics               numerics,
            java.io.BufferedWriter writer
            ) throws IOException {
        writer.write("// ==============================================\n");
        writer.write("// Please do not change the implementations below\n");
        writer.write("// =============================---==============\n");
        writer.write("#include \"" + solverName + ".h\"\n");
        writer.write( "#include \"kernels/aderdg/optimised/Kernels.h\"\n");
        writer.write("\n\n\n");
        
        switch(numerics) {
            case NONLINEAR:
                writer.write( "void " + _projectName + "::" + solverName + "::spaceTimePredictor( double* lQi, double* lFi, double* lQhi, double* lFhi, double* lQhbnd, double* lFhbnd, const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx, const double dt ) {\n");
                writer.write("   kernels::aderdg::optimised::picardLoop<flux>( lQi, lFi, luh, dx, dt );\n");
                writer.write("   kernels::aderdg::optimised::predictor( lQhi, lFhi, lQi, lFi );\n");
                writer.write("   kernels::aderdg::optimised::extrapolator( lQhbnd, lFhbnd, lQhi, lFhi );\n");
                writer.write("}\n");            
                break;
            case LINEAR:
                // TODO
                // Cauchy-Kowalewski
                break;
            default:
                break;
        }
        
        writer.write("\n\n\n");
        writer.write( "void " + _projectName + "::" + solverName + "::solutionUpdate(double* luh, const double* const lduh, const double dt) {\n");
        writer.write("   kernels::aderdg::optimised::solutionUpdate( luh, lduh, dt );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write( "void " + _projectName + "::" + solverName + "::volumeIntegral(double* lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
        writer.write("   kernels::aderdg::optimised::volumeIntegral( lduh, lFhi, dx );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write( "void " + _projectName + "::" + solverName + "::surfaceIntegral(double* lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
        writer.write("   kernels::aderdg::optimised::surfaceIntegral( lduh, lFhbnd, dx );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write( "void " + _projectName + "::" + solverName + "::riemannSolver(double* FL, double* FR, const double* const QL, const double* const QR, const double dt, const int normalNonZeroIndex) {\n");
        writer.write("   kernels::aderdg::optimised::riemannSolver<eigenvalues>( FL, FR, QL, QR, dt, normalNonZeroIndex );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write( "double " + _projectName + "::" + solverName + "::stableTimeStepSize(const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
        writer.write("   return kernels::aderdg::optimised::stableTimeStepSize<eigenvalues>( luh, dx );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write( "void " + _projectName + "::" + solverName + "::solutionAdjustment(double*  luh, const tarch::la::Vector<DIMENSIONS,double>&   center, const tarch::la::Vector<DIMENSIONS,double>&   dx, double  t, double  dt) {\n");
        writer.write("   kernels::aderdg::generic::solutionAdjustment<adjustedSolutionValues>( luh, center, dx, t, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
    }

    // @warning unused (please keep the code)
    private void writeADERDGSolverGeneratedImplementationForKernelEuler2d(
            String                 solverName,
            String                 numberOfVariables,
            String                 order,
            java.io.BufferedWriter writer
            ) throws IOException {
        writer.write("// ==============================================\n");
        writer.write("// Please do not change the implementations below\n");
        writer.write("// =============================---==============\n");
        writer.write("#include \"" + solverName + ".h\"\n");
        writer.write( "#include \"kernels/aderdg/optimised/Kernels.h\"\n");
        writer.write("\n\n\n");
        writer.write( "void " + _projectName + "::" + solverName + "::spaceTimePredictor( double* lQi, double* lFi, double* lQhi, double* lFhi, double* lQhbnd, double* lFhbnd, const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx, const double dt ) {\n");
        writer.write("   kernels::aderdg::optimised::picardLoop<flux>( lQi, lFi, luh, dx, dt );\n");
        writer.write("   kernels::aderdg::optimised::predictor( lQhi, lFhi, lQi, lFi );\n");
        writer.write("   kernels::aderdg::optimised::extrapolator( lQhbnd, lFhbnd, lQhi, lFhi );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write( "void " + _projectName + "::" + solverName + "::solutionUpdate(double* luh, const double* const lduh, const double dt) {\n");
        writer.write("   kernels::aderdg::optimised::solutionUpdate( luh, lduh, dt );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write( "void " + _projectName + "::" + solverName + "::volumeIntegral(double* lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
        writer.write("   kernels::aderdg::optimised::volumeIntegral( lduh, lFhi, dx );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write( "void " + _projectName + "::" + solverName + "::surfaceIntegral(double* lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
        writer.write("   kernels::aderdg::optimised::surfaceIntegral( lduh, lFhbnd, dx );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write( "void " + _projectName + "::" + solverName + "::riemannSolver(double* FL, double* FR, const double* const QL, const double* const QR, const double dt, const int normalNonZeroIndex) {\n");
        writer.write("   kernels::aderdg::optimised::riemannSolver<eigenvalues>( FL, FR, QL, QR, dt, normalNonZeroIndex );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write( "double " + _projectName + "::" + solverName + "::stableTimeStepSize(const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
        writer.write("   return kernels::aderdg::optimised::stableTimeStepSize<eigenvalues>( luh, dx );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write( "void " + _projectName + "::" + solverName + "::solutionAdjustment(double*  luh, const tarch::la::Vector<DIMENSIONS,double>&   center, const tarch::la::Vector<DIMENSIONS,double>&   dx, double  t, double  dt) {\n");
        writer.write("   kernels::aderdg::generic::solutionAdjustment<adjustedSolutionValues>( luh, center, dx, t, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
    }

*/

    private void invokeCodeGenerator(
            String                 solverName,
            String                 numberOfVariables,
            String                 order,
            Numerics               numerics
            ) throws IOException {
        
        String currentDirectory = System.getProperty("user.dir");
        java.nio.file.Path pathToCodeGenerator = java.nio.file.Paths.get(currentDirectory+"/Miscellaneous/CodeGenerator/Driver.py");
		    if(java.nio.file.Files.notExists(pathToCodeGenerator)) {
			    System.err.println("ERROR: Code generator not found. Can't generated optimised kernels.");
			    return;
		    }

		    String numericsParameter = "";
        switch(numerics) {
            case NONLINEAR:
                numericsParameter = "nonlinear";
                break;
            case LINEAR:
                numericsParameter = "linear";
                break;
            default:
                break;
        }		    

		    // set up the command to execute the code generator
		    String args          = " " + solverName                    + " "
		                               + numberOfVariables             + " "
		                               + order                         + " "
		                               + Integer.toString(_dimensions) + " "
		                               + numericsParameter             + " "
		                               + _microarchitecture            + " "
		                               + _pathToLibxsmm                + " "
		                               + "--precision=DP";  //double precision
		                               
		    String bashCommand   = "python " + pathToCodeGenerator + args ;

		    Runtime runtime = Runtime.getRuntime();

	      // execute the command line program
	      Process codeGenerator = runtime.exec(bashCommand);

	      // capture any output that is produced by the code generator and print it line-by-line
	      java.io.BufferedReader codeGeneratorsOutputReader = new java.io.BufferedReader(new java.io.InputStreamReader(codeGenerator.getInputStream()));
	      String line = "";
	      while((line = codeGeneratorsOutputReader.readLine()) != null) {
		        System.out.println(line);
	      }
    }
} 
