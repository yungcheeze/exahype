package eu.exahype;

import java.io.IOException;

import eu.exahype.analysis.DepthFirstAdapter;
import eu.exahype.node.AAderdgSolver;
import eu.exahype.node.AComputationalDomain;
import eu.exahype.node.AFiniteVolumesSolver;
import eu.exahype.node.ALimitingAderdgSolver;
import eu.exahype.node.AProfiling;
import eu.exahype.node.AProject;
import eu.exahype.node.PSolver;
import eu.exahype.solvers.Solver;
import eu.exahype.solvers.SolverFactory;

public class CreateSolverClasses extends DepthFirstAdapter {
  public Boolean valid = true;

  private DirectoryAndPathChecker _directoryAndPathChecker;

  private String _projectName;

  private String _microarchitecture;

  private java.util.List<String> _supportedMicroarchitectures;

  private java.util.Set<String>  _definedSolvers;

  private String _pathToLibxsmm;

  private int _dimensions;

  private boolean _enableProfiler;

  public CreateSolverClasses(DirectoryAndPathChecker directoryAndPathChecker) {
    _directoryAndPathChecker = directoryAndPathChecker;
    _supportedMicroarchitectures =
        java.util.Arrays.asList("wsm", "snb", "hsw", "knc", "knl", "noarch");
    _enableProfiler = false;
  }

  @Override
  public void inAProject(AProject node) {
    _projectName     = node.getName().toString().trim();
    _definedSolvers  = new java.util.HashSet<String>();

    if (node.getSolver().size() == 0) {
      System.out.println("there are no solvers in the specification file ... nothing to be done");
    }

    // Only one optimised solver can be used (optimised kernel would be overwritten by the latest solver otherwise)
    if (node.getSolver().size() > 1) {
      int optimisedCount = 0;
      for(PSolver psolver : node.getSolver()) {
        if(psolver instanceof AAderdgSolver) {
          AAderdgSolver asolver = (AAderdgSolver) psolver;
          if(    asolver.getKernel().toString().trim().equals( eu.exahype.solvers.OptimisedFluxesNonlinearADER_DGinC.Identifier )
              || asolver.getKernel().toString().trim().equals( eu.exahype.solvers.OptimisedFluxesLinearADER_DGinC.Identifier )
              ){
            optimisedCount++;
          }
        }
      }
      if(optimisedCount > 1) {
        System.err.println("ERROR: Only one optimised solver can be used at a time. Currently "+optimisedCount+" are defined.");
        valid = false;
        return;
      }
    }

    _microarchitecture = node.getArchitecture().toString().trim().toLowerCase();
    if (!_supportedMicroarchitectures.contains(_microarchitecture)) {
      System.out.println("Unknown architecture specified ... fallback solution \"noarch\" taken");
      _microarchitecture = "noarch";
    }
  }

  @Override
  public void inAPaths(eu.exahype.node.APaths node) {
    if (node.getLibxsmmPath() == null) {
      // attribute 'libxsmm-path' did not occur in spec file
      _pathToLibxsmm = "";
    } else {
      _pathToLibxsmm = node.getLibxsmmPath().toString().trim();
    }
  };

  @Override
  public void inAComputationalDomain(AComputationalDomain node) {
    _dimensions = Integer.parseInt( node.getDimension().toString().trim() );
    if (_dimensions!=2 && _dimensions!=3) {
      System.err.println( "ERROR: dimension has to be either 2 or 3.");
    }
  }


  @Override
  public void inAProfiling(AProfiling node) {
    _enableProfiler = !node.getProfiler().toString().trim().equals("NoOpProfiler");
  };

  @Override
  public void inAAderdgSolver(AAderdgSolver node) {
    String solverName = node.getName().toString().trim();

    if (_definedSolvers.contains(solverName)) {
      System.err.println( "ERROR: Solver " + solverName + " multiply defined" );
      valid = false;
    }
    else {
      _definedSolvers.add(solverName);
    }

    java.io.File userPDEFile = null;
    java.io.File userTypesDefFile = null;

    boolean isFortran = false;
    if (node.getLanguage().getText().trim().equals("C")) {
      isFortran = false;
    } else if (node.getLanguage().getText().trim().equals("Fortran")) {
      isFortran = true;
      userPDEFile =
          new java.io.File(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/PDE.f90");
      userTypesDefFile = new java.io.File(
          _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/typesDef.f90");
    } else {
      System.err.println("ERROR: unknown language for solver " + node.getName().getText()
          + ". Supported languages are C and Fortran");
      valid = false;
      return;
    }

    String  kernel             = node.getKernel().toString().trim();
    int     numberOfVariables  = Integer.parseInt(node.getVariables().toString().trim());
    int     numberOfParameters = Integer.parseInt(node.getParameters().toString().trim());
    int     order              = Integer.parseInt(node.getOrder().toString().trim());
    boolean hasConstants       = node.getConstants()!=null;

    if (numberOfParameters != 0) {
      System.err.println("ERROR: At the moment, parameters are not yet supported. " + 
          " Please add the parameters as additional quantities to your PDE formulation.");
      valid = false;
      return;
    }

    if (order < 1 || order > 9) {
      System.err.println("ERROR: Only polynomial degrees of 1..9 are supported.");
      valid = false;
      return;
    }
    
    SolverFactory solverFactory = new SolverFactory(_dimensions, _enableProfiler, _microarchitecture, _pathToLibxsmm);
    eu.exahype.solvers.Solver solver = solverFactory.createADERDGSolver(
        kernel, isFortran, numberOfVariables, numberOfParameters, order, hasConstants);

    if (solver == null) {
      System.err.println("creation solver " + solverName + " ... failed as kernel " + kernel
          + " for language " + node.getLanguage().getText().trim() + " is not supported");
      valid = false;
      return;
    }

    //
    // Write the files
    //
    try {
      tryWriteSolverHeader(solver, solverName);
      tryWriteSolverUserImplementation(solver, solverName);
      tryWriteSolverGeneratedImplementation(solver, solverName);
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      valid = false;
    }
  }

  @Override
  public void inAFiniteVolumesSolver(AFiniteVolumesSolver node) {
    String solverName = node.getName().toString().trim();

    if (_definedSolvers.contains(solverName)) {
      System.err.println( "ERROR: Solver " + solverName + " multiply defined" );
      valid = false;
    }
    else {
      _definedSolvers.add(solverName);
    }

    java.io.File userPDEFile = null;
    java.io.File userTypesDefFile = null;


    boolean isFortran = false;
    if (node.getLanguage().getText().trim().equals("C")) {
      isFortran = false;
    } else if (node.getLanguage().getText().trim().equals("Fortran")) {
      isFortran = true;
      userPDEFile =
          new java.io.File(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/PDE.f90");
      userTypesDefFile = new java.io.File(
          _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/typesDef.f90");
    } else {
      System.err.println("ERROR: unknown language for solver " + node.getName().getText()
          + ". Supported languages are C and Fortran");
      valid = false;
      return;
    }

    String kernel = node.getKernel().toString().trim();

    int numberOfVariables  = Integer.parseInt(node.getVariables().toString().trim());
    int numberOfParameters = Integer.parseInt(node.getParameters().toString().trim());
    int patchSize          = Integer.parseInt(node.getPatchSize().toString().trim());
    boolean hasConstants   = node.getConstants()!=null;

    if (numberOfParameters != 0) {
      System.err.println("ERROR: At the moment, parameters are not yet supported. " + 
          " Please add the parameters as additional quantities to your PDE formulation.");
      valid = false;
      return;
    }
    
    SolverFactory solverFactory = new SolverFactory(_dimensions, _enableProfiler, _microarchitecture, _pathToLibxsmm);
    Solver solver = solverFactory.createFiniteVolumesSolver(
        kernel,isFortran,numberOfVariables,numberOfParameters,patchSize,hasConstants);

    if (solver == null) {
      System.err.println("creation solver " + solverName + " ... failed as kernel " + kernel
          + " for language " + node.getLanguage().getText().trim() + " is not supported");
      valid = false;
      return;
    }

    //
    // Write the files
    //
    try {
      tryWriteSolverHeader(solver, solverName);
      
      tryWriteSolverUserImplementation(solver, solverName);
      
      tryWriteSolverGeneratedImplementation(solver, solverName);
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      valid = false;
    }
  }

  @Override
  public void inALimitingAderdgSolver(ALimitingAderdgSolver node) {
    String solverName = node.getName().toString().trim();

    if (_definedSolvers.contains(solverName)) {
      System.err.println( "ERROR: Solver " + solverName + " multiply defined" );
      valid = false;
    }
    else {
      _definedSolvers.add(solverName);
    }

    java.io.File userPDEFile = null; // TODO(Dominic): Fortran specifics; not used yet
    java.io.File userTypesDefFile = null;
    

    boolean isFortran = false;
    if (node.getLanguage().getText().trim().equals("C")) {
      isFortran = false;
    } else if (node.getLanguage().getText().trim().equals("Fortran")) {
      isFortran = true;
      userPDEFile =
          new java.io.File(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/PDE.f90");
      userTypesDefFile = new java.io.File(
          _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/typesDef.f90");
    } else {
      System.err.println("ERROR: unknown language for solver " + node.getName().getText()
          + ". Supported languages are C and Fortran");
      valid = false;
      return;
    }

    String  kernel             = node.getKernel().toString().trim();
    int     numberOfVariables  = Integer.parseInt(node.getVariables().toString().trim());
    int     numberOfParameters = Integer.parseInt(node.getParameters().toString().trim());
    int     order              = Integer.parseInt(node.getOrder().toString().trim());
    int     patchSize          = 2*order+1;
    boolean hasConstants       = node.getConstants()!=null;
    
    String  limiterKernel      = node.getKernelLimiter().toString().trim();

    if (numberOfParameters != 0) {
      System.err.println("ERROR: At the moment, parameters are not yet supported. " + 
          " Please add the parameters as additional quantities to your PDE formulation.");
      valid = false;
      return;
    }
    
    SolverFactory solverFactory = new SolverFactory(_dimensions, _enableProfiler, _microarchitecture, _pathToLibxsmm);
    Solver solver  = solverFactory.createADERDGSolver(
        kernel,isFortran,numberOfVariables,numberOfParameters,order,hasConstants);
    Solver limiter = solverFactory.createFiniteVolumesSolver(
        limiterKernel,isFortran,numberOfVariables,numberOfParameters,patchSize,hasConstants);

    if (solver == null || limiter == null) {
      System.err.println("creation solver " + solverName + " ... failed as kernel " + kernel
          + " for language " + node.getLanguage().getText().trim() + " is not supported");
      valid = false;
      return;
    }

    //
    // Write the files
    //
    try {
      tryWriteSolverHeader(solver, solverName+"_ADERDG");
      tryWriteSolverHeader(limiter, solverName+"_FV");

      tryWriteSolverUserImplementation(solver,solverName+"_ADERDG");
      tryWriteSolverUserImplementation(limiter,solverName+"_FV");

      tryWriteSolverGeneratedImplementation(solver,solverName+"_ADERDG");
      tryWriteSolverGeneratedImplementation(limiter,solverName+"_FV");
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      valid = false;
    }
  }
  
  private void tryWriteSolverHeader(Solver solver,String solverName) throws IOException {
    java.io.File solverHeaderFile = new java.io.File(
        _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/" + solverName + ".h");
    
    if (solverHeaderFile.exists()) {
      System.out.println("create header of solver " + solverName + " ... header "
          + solverHeaderFile.getAbsoluteFile()
          + " does exist already. Remove to allow toolkit to regenerate it (changes will be lost)");
    } else {
      java.io.BufferedWriter headerWriter =
          new java.io.BufferedWriter(new java.io.FileWriter(solverHeaderFile));
      solver.writeHeader(headerWriter, solverName, _projectName);
      System.out.println("create header of solver " + solverName + " ... ok");
      headerWriter.close();
    }
  }
  
  private void tryWriteSolverUserImplementation(Solver solver, String solverName) throws IOException {
    java.io.File solverUserImplementationFile = new java.io.File(
        _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/" + solverName + ".cpp");
    
    if (solverUserImplementationFile.exists()) {
      System.out.println("user's implementation file of solver " + solverName
          + " ... does exist already. Is not overwritten");
    } else {
      java.io.BufferedWriter userImplementationWriter =
          new java.io.BufferedWriter(new java.io.FileWriter(solverUserImplementationFile));
      solver.writeUserImplementation(userImplementationWriter, solverName, _projectName);
      System.out.println(
          "create user implementation template of solver " + solverName + " ... please complete");
      userImplementationWriter.close();
    }
  }
  
  private void tryWriteSolverGeneratedImplementation(Solver solver, String solverName) throws IOException {
    java.io.File solverGeneratedImplementationFile = new java.io.File(
        _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/" + solverName + "_generated.cpp");
    
    if (solverGeneratedImplementationFile.exists()) {
      System.out.println("generated implementation file of solver " + solverName
          + " ... does exist already. Is overwritten");
    }

    java.io.BufferedWriter generatedImplementationWriter =
        new java.io.BufferedWriter(new java.io.FileWriter(solverGeneratedImplementationFile));
    solver.writeGeneratedImplementation(generatedImplementationWriter, solverName, _projectName);
    System.out.println("create generated implementation of solver " + solverName + " ... ok");
    generatedImplementationWriter.close();
  }
}
