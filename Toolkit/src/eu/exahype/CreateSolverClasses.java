package eu.exahype;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

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
import eu.exahype.FileSearch;

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
  
  public static Map<String,Integer> getVariables(PSolver node) {
    String variablesAsString = null;
    
    if (node instanceof AAderdgSolver) {
      variablesAsString = ((AAderdgSolver) node).getVariables().getText();
    } else if (node instanceof ALimitingAderdgSolver) {
      variablesAsString = ((ALimitingAderdgSolver) node).getVariables().getText();
    } else if (node instanceof AFiniteVolumesSolver) {
      variablesAsString = ((AFiniteVolumesSolver) node).getVariables().getText();
    } 
//    else if (node instanceof AAderdgWithVariablesListSolver) {
//      variablesAsString = ((AAderdgWithVariablesListSolver) node).getVariables().getText();
//    } else if (node instanceof ALimitingAderdgWithVariablesListSolver) {
//      variablesAsString = ((ALimitingAderdgWithVariablesListSolver) node).getVariables().getText();
//    } else if (node instanceof AFiniteVolumesWithVariablesListSolver) {
//      variablesAsString = ((AFiniteVolumesWithVariablesListSolver) node).getVariables().getText();
//    } 
    else {
      System.out.println("ERROR: I do not know how to handle solver type "+node.getClass().toString()+"!");
      System.exit(1);
    }
    
    try { // the user only gave us a number, e.g., 5, instead of a list, e.g., v0:1,v1:3,v2:3.
      int numberOfVariables = Integer.parseInt(variablesAsString);
      
      Map<String,Integer> map = new HashMap<String, Integer>(1);
      map.put("Q", numberOfVariables);
      return map;
    } catch (NumberFormatException exception) { // the user gave us a list       
      String[] variables = variablesAsString.split(",");
      Map<String,Integer> map = new HashMap<String, Integer>(variables.length);
      
      for (String variable : variables) {
        String[] identifierAndQuantity = variable.split(":");
        
        String identifier = identifierAndQuantity[0].trim();
        // TODO(Dominic): Remove from here and move to Acess object generation code
        //        { 
//          identifier      = identifier.substring(0,1).toUpperCase() + identifier.substring(1).toLowerCase();
//        } 
        try {
          int dimension    = Integer.parseInt(identifierAndQuantity[1].trim());
          
          if (dimension <= 0) {
            System.out.println("ERROR: Quantity specifier of '"+identifier+"' is not a positive integer!");
            System.exit(1);
          }
          
          map.put(identifier,dimension);
          
          System.out.println("Found variable "+identifier+" with "+dimension+" elements."); // TODO(Dominic): Comment in for debugging purposes
          
        } catch (NumberFormatException exception2) { 
          System.out.println("ERROR: Quantity specifier of '"+identifier+"' is not a positive integer!");
          System.exit(1);
        }
      }
      return map;
    }
  }
  
  public static int getNumberOfVariables(Map<String,Integer> variables) {
    int numberOfVariables = 0;
    for (String key : variables.keySet()) {
      numberOfVariables += variables.get(key);
    }
    
    System.out.println("Total number of state variables is "+numberOfVariables); // TODO(Dominic): Comment in for debugging purposes
    
    return numberOfVariables;
  }

  public CreateSolverClasses(DirectoryAndPathChecker directoryAndPathChecker) {
    _directoryAndPathChecker = directoryAndPathChecker;
    _supportedMicroarchitectures =
        java.util.Arrays.asList("wsm", "snb", "hsw", "knc", "knl", "noarch");
    _enableProfiler = false;
  }

  @Override
  public void inAProject(AProject node) {
    _projectName     = node.getName().getText();
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
          if(    asolver.getKernel().getText().equals( eu.exahype.solvers.OptimisedFluxesNonlinearADER_DGinC.Identifier )
              || asolver.getKernel().getText().equals( eu.exahype.solvers.OptimisedFluxesLinearADER_DGinC.Identifier )
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

    _microarchitecture = node.getArchitecture().getText().toLowerCase();
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
      _pathToLibxsmm = node.getLibxsmmPath().getText();
    }
  };

  @Override
  public void inAComputationalDomain(AComputationalDomain node) {
    _dimensions = Integer.parseInt( node.getDimension().getText() );
    if (_dimensions!=2 && _dimensions!=3) {
      System.err.println( "ERROR: dimension has to be either 2 or 3.");
    }
  }


  @Override
  public void inAProfiling(AProfiling node) {
    _enableProfiler = !node.getProfiler().getText().equals("NoOpProfiler");
  };

  @Override
  public void inAAderdgSolver(AAderdgSolver node) {
    String solverName = node.getName().getText();

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
          FileSearch.relocatableFile(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/PDE.f90");
      userTypesDefFile = FileSearch.relocatableFile(
          _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/typesDef.f90");
    } else {
      System.err.println("ERROR: unknown language for solver " + node.getName().getText()
          + ". Supported languages are C and Fortran");
      valid = false;
      return;
    }

    String  kernel                = node.getKernel().getText();
    Map<String,Integer> variables = getVariables(node);
    int     numberOfVariables     = getNumberOfVariables(variables);
    int     numberOfParameters    = Integer.parseInt(node.getParameters().getText());
    int     order                 = Integer.parseInt(node.getOrder().getText());
    boolean hasConstants          = node.getConstants()!=null;

    if (numberOfParameters != 0) {
      System.err.println("ERROR: At the moment, parameters are not supported. " + 
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
    String solverName = node.getName().getText();

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
          FileSearch.relocatableFile(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/PDE.f90");
      userTypesDefFile = FileSearch.relocatableFile(
          _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/typesDef.f90");
    } else {
      System.err.println("ERROR: unknown language for solver " + node.getName().getText()
          + ". Supported languages are C and Fortran");
      valid = false;
      return;
    }

    String kernel = node.getKernel().getText();

    Map<String,Integer> variables = getVariables(node);
    int     numberOfVariables     = getNumberOfVariables(variables);
    int numberOfParameters        = Integer.parseInt(node.getParameters().getText());
    int patchSize                 = Integer.parseInt(node.getPatchSize().getText());
    boolean hasConstants          = node.getConstants()!=null;

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
    String solverName = node.getName().getText();

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
          FileSearch.relocatableFile(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/PDE.f90");
      userTypesDefFile = FileSearch.relocatableFile(
          _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/typesDef.f90");
    } else {
      System.err.println("ERROR: unknown language for solver " + node.getName().getText()
          + ". Supported languages are C and Fortran");
      valid = false;
      return;
    }

    String  kernel                = node.getKernel().getText();
    Map<String,Integer> variables = getVariables(node);
    int     numberOfVariables     = getNumberOfVariables(variables);
    int     numberOfParameters    = Integer.parseInt(node.getParameters().getText());
    int     order                 = Integer.parseInt(node.getOrder().getText());
    int     patchSize             = 2*order+1;
    boolean hasConstants          = node.getConstants()!=null;
    
    String  limiterKernel      = node.getKernelLimiter().getText();

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
  
//  /**
//   * The function behaviour is identical to inAAderdgSolver.
//   */
//  @Override
//  public void inAAderdgWithVariablesListSolver(AAderdgWithVariablesListSolver node) {
//    String solverName = node.getName().getText();
//
//    if (_definedSolvers.contains(solverName)) {
//      System.err.println( "ERROR: Solver " + solverName + " multiply defined" );
//      valid = false;
//    }
//    else {
//      _definedSolvers.add(solverName);
//    }
//
//    java.io.File userPDEFile = null;
//    java.io.File userTypesDefFile = null;
//
//    boolean isFortran = false;
//    if (node.getLanguage().getText().trim().equals("C")) {
//      isFortran = false;
//    } else if (node.getLanguage().getText().trim().equals("Fortran")) {
//      isFortran = true;
//      userPDEFile =
//          FileSearch.relocatableFile(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/PDE.f90");
//      userTypesDefFile = FileSearch.relocatableFile(
//          _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/typesDef.f90");
//    } else {
//      System.err.println("ERROR: unknown language for solver " + node.getName().getText()
//          + ". Supported languages are C and Fortran");
//      valid = false;
//      return;
//    }
//
//    String  kernel                = node.getKernel().getText();
//    Map<String,Integer> variables = getVariables(node);
//    int     numberOfVariables     = getNumberOfVariables(variables);
//    int     numberOfParameters    = Integer.parseInt(node.getParameters().getText());
//    int     order                 = Integer.parseInt(node.getOrder().getText());
//    boolean hasConstants          = node.getConstants()!=null;
//
//    if (numberOfParameters != 0) {
//      System.err.println("ERROR: At the moment, parameters are not supported. " + 
//          " Please add the parameters as additional quantities to your PDE formulation.");
//      valid = false;
//      return;
//    }
//
//    if (order < 1 || order > 9) {
//      System.err.println("ERROR: Only polynomial degrees of 1..9 are supported.");
//      valid = false;
//      return;
//    }
//    
//    SolverFactory solverFactory = new SolverFactory(_dimensions, _enableProfiler, _microarchitecture, _pathToLibxsmm);
//    eu.exahype.solvers.Solver solver = solverFactory.createADERDGSolver(
//        kernel, isFortran, numberOfVariables, numberOfParameters, order, hasConstants);
//
//    if (solver == null) {
//      System.err.println("creation solver " + solverName + " ... failed as kernel " + kernel
//          + " for language " + node.getLanguage().getText().trim() + " is not supported");
//      valid = false;
//      return;
//    }
//
//    //
//    // Write the files
//    //
//    try {
//      tryWriteSolverHeader(solver, solverName);
//      tryWriteSolverUserImplementation(solver, solverName);
//      tryWriteSolverGeneratedImplementation(solver, solverName);
//    } catch (Exception exc) {
//      System.err.println("ERROR: " + exc.toString());
//      valid = false;
//    }
//  }
//  
//  /**
//   * The function behaviour is identical to inAFiniteVolumesSolver.
//   */
//  @Override
//  public void inAFiniteVolumesWithVariablesListSolver(AFiniteVolumesWithVariablesListSolver node) {
//    String solverName = node.getName().getText();
//
//    if (_definedSolvers.contains(solverName)) {
//      System.err.println( "ERROR: Solver " + solverName + " multiply defined" );
//      valid = false;
//    }
//    else {
//      _definedSolvers.add(solverName);
//    }
//
//    java.io.File userPDEFile = null;
//    java.io.File userTypesDefFile = null;
//
//
//    boolean isFortran = false;
//    if (node.getLanguage().getText().trim().equals("C")) {
//      isFortran = false;
//    } else if (node.getLanguage().getText().trim().equals("Fortran")) {
//      isFortran = true;
//      userPDEFile =
//          FileSearch.relocatableFile(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/PDE.f90");
//      userTypesDefFile = FileSearch.relocatableFile(
//          _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/typesDef.f90");
//    } else {
//      System.err.println("ERROR: unknown language for solver " + node.getName().getText()
//          + ". Supported languages are C and Fortran");
//      valid = false;
//      return;
//    }
//
//    String kernel = node.getKernel().getText();
//
//    Map<String,Integer> variables = getVariables(node);
//    int     numberOfVariables     = getNumberOfVariables(variables);
//    int numberOfParameters        = Integer.parseInt(node.getParameters().getText());
//    int patchSize                 = Integer.parseInt(node.getPatchSize().getText());
//    boolean hasConstants          = node.getConstants()!=null;
//
//    if (numberOfParameters != 0) {
//      System.err.println("ERROR: At the moment, parameters are not yet supported. " + 
//          " Please add the parameters as additional quantities to your PDE formulation.");
//      valid = false;
//      return;
//    }
//    
//    SolverFactory solverFactory = new SolverFactory(_dimensions, _enableProfiler, _microarchitecture, _pathToLibxsmm);
//    Solver solver = solverFactory.createFiniteVolumesSolver(
//        kernel,isFortran,numberOfVariables,numberOfParameters,patchSize,hasConstants);
//
//    if (solver == null) {
//      System.err.println("creation solver " + solverName + " ... failed as kernel " + kernel
//          + " for language " + node.getLanguage().getText().trim() + " is not supported");
//      valid = false;
//      return;
//    }
//
//    //
//    // Write the files
//    //
//    try {
//      tryWriteSolverHeader(solver, solverName);
//      
//      tryWriteSolverUserImplementation(solver, solverName);
//      
//      tryWriteSolverGeneratedImplementation(solver, solverName);
//    } catch (Exception exc) {
//      System.err.println("ERROR: " + exc.toString());
//      valid = false;
//    }
//  }
//  
//  /**
//   * The function behaviour is identical to inALimitingAderdgWithVariablesListSolver.
//   */
//  @Override
//  public void inALimitingAderdgWithVariablesListSolver(ALimitingAderdgWithVariablesListSolver node) {
//    String solverName = node.getName().getText();
//
//    if (_definedSolvers.contains(solverName)) {
//      System.err.println( "ERROR: Solver " + solverName + " multiply defined" );
//      valid = false;
//    }
//    else {
//      _definedSolvers.add(solverName);
//    }
//
//    java.io.File userPDEFile = null; // TODO(Dominic): Fortran specifics; not used yet
//    java.io.File userTypesDefFile = null;
//    
//
//    boolean isFortran = false;
//    if (node.getLanguage().getText().trim().equals("C")) {
//      isFortran = false;
//    } else if (node.getLanguage().getText().trim().equals("Fortran")) {
//      isFortran = true;
//      userPDEFile =
//          FileSearch.relocatableFile(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/PDE.f90");
//      userTypesDefFile = FileSearch.relocatableFile(
//          _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/typesDef.f90");
//    } else {
//      System.err.println("ERROR: unknown language for solver " + node.getName().getText()
//          + ". Supported languages are C and Fortran");
//      valid = false;
//      return;
//    }
//
//    String  kernel                = node.getKernel().getText();
//    Map<String,Integer> variables = getVariables(node);
//    int     numberOfVariables     = getNumberOfVariables(variables);
//    int     numberOfParameters    = Integer.parseInt(node.getParameters().getText());
//    int     order                 = Integer.parseInt(node.getOrder().getText());
//    int     patchSize             = 2*order+1;
//    boolean hasConstants          = node.getConstants()!=null;
//    
//    String  limiterKernel      = node.getKernelLimiter().getText();
//
//    if (numberOfParameters != 0) {
//      System.err.println("ERROR: At the moment, parameters are not yet supported. " + 
//          " Please add the parameters as additional quantities to your PDE formulation.");
//      valid = false;
//      return;
//    }
//    
//    SolverFactory solverFactory = new SolverFactory(_dimensions, _enableProfiler, _microarchitecture, _pathToLibxsmm);
//    Solver solver  = solverFactory.createADERDGSolver(
//        kernel,isFortran,numberOfVariables,numberOfParameters,order,hasConstants);
//    Solver limiter = solverFactory.createFiniteVolumesSolver(
//        limiterKernel,isFortran,numberOfVariables,numberOfParameters,patchSize,hasConstants);
//
//    if (solver == null || limiter == null) {
//      System.err.println("creation solver " + solverName + " ... failed as kernel " + kernel
//          + " for language " + node.getLanguage().getText().trim() + " is not supported");
//      valid = false;
//      return;
//    }
//
//    //
//    // Write the files
//    //
//    try {
//      tryWriteSolverHeader(solver, solverName+"_ADERDG");
//      tryWriteSolverHeader(limiter, solverName+"_FV");
//
//      tryWriteSolverUserImplementation(solver,solverName+"_ADERDG");
//      tryWriteSolverUserImplementation(limiter,solverName+"_FV");
//
//      tryWriteSolverGeneratedImplementation(solver,solverName+"_ADERDG");
//      tryWriteSolverGeneratedImplementation(limiter,solverName+"_FV");
//    } catch (Exception exc) {
//      System.err.println("ERROR: " + exc.toString());
//      valid = false;
//    }
//  }
  
  private void tryWriteSolverHeader(Solver solver,String solverName) throws IOException {
    java.io.File solverHeaderFile = FileSearch.relocatableFile(
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
    java.io.File solverUserImplementationFile = FileSearch.relocatableFile(
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
    java.io.File solverGeneratedImplementationFile = FileSearch.relocatableFile(
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
