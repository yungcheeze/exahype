package eu.exahype;

import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;

public class CodeGeneratorHelper {
  
  //configuration parameters
  private static String OPT_KERNEL_PATH_PREFIX = "kernels";                  //starts from application root
  private static String CODEGENERATOR_PATH     = "CodeGenerator/Driver.py";  //starts from exahype root (ExaHyPE-Engine)
  
  //Singleton pattern to be able to access the instance in solvers
  private static volatile CodeGeneratorHelper instance = null;

  private CodeGeneratorHelper() {
    _optKernelPaths = new HashSet<String>();
  }

  public static CodeGeneratorHelper getInstance() {
      if (instance == null) {
          synchronized(CodeGeneratorHelper.class) {
              if (instance == null) {
                  instance = new CodeGeneratorHelper();
              }
          }
      }
      return instance;
  }
  
  
  //Internal states
  private Collection<String> _optKernelPaths;
  private String _pathToLibxsmm = null;
  private String _pathToApplication = null;
  
  //setter
  public void setPaths(DirectoryAndPathChecker directoryAndPathChecker) {
    try{
      _pathToLibxsmm = directoryAndPathChecker.libxsmmPath.getCanonicalPath();
      _pathToApplication = directoryAndPathChecker.outputDirectory.getPath();
    } catch(IOException e) {} 
  }
  
  //getter
  public Collection<String> getOptKernelPaths() {
    return _optKernelPaths;
  }
  
  public String invokeCodeGenerator(String projectName, String solverName, int numberOfUnknowns, int numberOfParameters, int order,
      boolean isLinear, int dimensions, String microarchitecture, boolean enableDeepProfiler, boolean useFlux, boolean useSource, boolean useNCP, boolean noTimeAveraging)
      throws IOException {
        
    if(_pathToLibxsmm == null) {
      System.err.println("ERROR: Path to Libxsmm for the Code generator not found");
      throw new IOException();
    }
    
    if(_pathToApplication == null) {
      System.err.println("ERROR: Path to the application for the Code generator not found");
      throw new IOException();
    }
        
        
    java.io.File pathToCodeGenerator_File =
        new java.io.File(CodeGeneratorHelper.CODEGENERATOR_PATH);
    if (!pathToCodeGenerator_File.exists()) {
      System.err.println("ERROR: Code generator not found. Can't generate optimised kernels. Path: " + pathToCodeGenerator_File.getCanonicalPath());
      throw new IOException();
    }
    String pathToCodeGenerator = pathToCodeGenerator_File.getCanonicalPath();
    
    String numericsParameter = isLinear ? "linear" : "nonlinear";
    String options = (enableDeepProfiler ? "--deepProfiling " : "") + (useFlux ? "--useFlux " : "") + (useSource ? "--useSource " : "") + (useNCP ? "--useNCP " : "") + (noTimeAveraging ? "--noTimeAveraging " : "");
    

    // set up the command to execute the code generator
    String args =   " " + projectName + "::" + solverName 
                  + " " + numberOfUnknowns 
                  + " " + order 
                  + " " + dimensions 
                  + " " + numericsParameter 
                  + " " + microarchitecture 
                  + " " + options; 

    String optKernelPath = (new java.io.File(OPT_KERNEL_PATH_PREFIX,solverName)).getPath();
    
    String bashCommand = "env python3 "  + pathToCodeGenerator + " " + _pathToLibxsmm + " " + _pathToApplication + " " + optKernelPath + args;

    Runtime runtime = Runtime.getRuntime();
    System.out.println("Codegenerator command line: "+bashCommand);
    // execute the command line program
    Process codeGenerator = runtime.exec(bashCommand);

    // capture any output that is produced by the code generator and print it line-by-line
    java.io.InputStream stdout = codeGenerator.getInputStream();
    java.io.BufferedReader stdoutReader =
        new java.io.BufferedReader(new java.io.InputStreamReader(stdout));
    String line = "";
    while ((line = stdoutReader.readLine()) != null) {
      System.out.println("CodeGenerator: " + line);
    }
    java.io.InputStream stderr = codeGenerator.getErrorStream();
    java.io.BufferedReader stderrReader =
        new java.io.BufferedReader(new java.io.InputStreamReader(stderr));
    while ((line = stderrReader.readLine()) != null) {
      System.out.println("CodeGenerator: " + line);
    }

    // in order to stop further toolkit execution if the code generator fails,
    // explicitly wait for the process
    try {
        int exitValue = codeGenerator.waitFor();
        if(exitValue != 0) {
            System.err.println("ERROR: Code Generator failed with exit value " + exitValue);
            throw new IOException(); // <- also done in line 186. This is abusing the exception system.
        }
    } catch(InterruptedException e) {
        System.err.println("This is very bad. I don't know what's going on.");
        throw new IOException();
    }
    
     _optKernelPaths.add(optKernelPath);
    
    return optKernelPath;
    
  } // invokeCodeGenerator
  
  public static String getOptKernelPath(String key) {
    return OPT_KERNEL_PATH_PREFIX+"/"+key;
  }

}
