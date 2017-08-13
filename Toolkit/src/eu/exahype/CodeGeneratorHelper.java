package eu.exahype;

import java.io.IOException;
import java.nio.file.*;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.Map;
import java.util.HashMap;
import java.util.Collection;
import java.util.HashSet;
import java.util.Collections;

public class CodeGeneratorHelper {
  
  //configuration parameters
  private static String INTERNAL_EXAHYPE_PATH  = "ExaHyPE";                  //starts from root (ExaHyPE-Engine)
  private static String OPT_KERNEL_PATH_PREFIX = "kernels/aderdg/optimised"; //starts from root + INTERNAL_EXAHYPE_PATH
  private static String CODEGENERATOR_PATH     = "CodeGenerator/Driver.py";  //starts from root (ExaHyPE-Engine)
  
  //Singleton pattern to be able to access the instance in solvers
  private static volatile CodeGeneratorHelper instance = null;

  private CodeGeneratorHelper() {
    _optDirectories = new HashMap<String, Collection<String>>();
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
  private Map<String, Collection<String>> _optDirectories;
  
  //getter
  public Collection<String> getOptKernelPaths(String projectName) {
    if(_optDirectories.containsKey(projectName) && _optDirectories.get(projectName) != null) {
      return _optDirectories.get(projectName);
    } else {
      return Collections.<String>emptySet(); //return an immutable empty Collection<String>
    }
  }
  
  private void addPathToMap(String projectName, String path) {
    if(_optDirectories.containsKey(projectName) && _optDirectories.get(projectName) != null) {
      _optDirectories.get(projectName).add(path);
    } else {
      HashSet<String> set = new HashSet<String>();
      set.add(path);
      _optDirectories.put(projectName, set);
    }
  }
  
  
  public String invokeCodeGenerator(String projectName, String solverFullName, int numberOfUnknowns, int numberOfParameters, int order,
      boolean isLinear, int dimensions, String microarchitecture, String pathToLibxsmm, boolean enableDeepProfiler, boolean useFlux, boolean useSource, boolean useNCP, boolean noTimeAveraging)
      throws IOException {
    String currentDirectory = System.getProperty("user.dir");
    java.io.File pathToCodeGenerator =
        new java.io.File(currentDirectory + '/' + CodeGeneratorHelper.CODEGENERATOR_PATH);
    if (!pathToCodeGenerator.exists()) {
      System.err.println("ERROR: Code generator not found. Can't generate optimised kernels. Path: " + pathToCodeGenerator.toString());
      throw new IOException();
    }
    
    if(pathToLibxsmm == null || pathToLibxsmm.isEmpty()) {
      System.err.println("ERROR: Libxsmm path not specified");
      throw new IOException();
    }
    
/*    java.io.File pathToLibxsmmMakefile = //To test if the libxsmm folder is correct
        new java.io.File(java.nio.file.Paths.get(currentDirectory,pathToLibxsmm,"Makefile").toString());
    if (!pathToLibxsmmMakefile.exists()) {
      System.err.println("ERROR: Libxsmm makefile not found. Can't generate optimised kernels. Path: " + pathToLibxsmmMakefile.toString());
      throw new IOException();
    }
*/
    String numericsParameter = isLinear ? "linear" : "nonlinear";
    String options = (enableDeepProfiler ? "--deepProfiling " : "") + (useFlux ? "--useFlux " : "") + (useSource ? "--useSource " : "") + (useNCP ? "--useNCP " : "") + (noTimeAveraging ? "--noTimeAveraging " : "");
    

    // set up the command to execute the code generator
    String args = " " + solverFullName + " " + numberOfUnknowns + " " + order + " "
        + dimensions + " " + numericsParameter + " " + microarchitecture + " "
        + currentDirectory + "/"  + pathToLibxsmm + " " + options; 

    String optKernelPath = getOptKernelPath(args);
        
    String bashCommand = "env python3 "  + pathToCodeGenerator  + " " + optKernelPath + args;

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
    
    addPathToMap(projectName, optKernelPath);
    
    return optKernelPath;
    
  } // invokeCodeGenerator
  
  public static String getOptKernelPath(String key) {
    //TODO JMG
    
    return OPT_KERNEL_PATH_PREFIX+"/test_TODOJMG";
  }
  
  //remove all the generated opt. kernels
  public static void cleanAll() throws IOException {
    System.out.println("Cleaning directory "+Paths.get(INTERNAL_EXAHYPE_PATH, OPT_KERNEL_PATH_PREFIX).toString());
    cleanDirectory(Paths.get(System.getProperty("user.dir"), INTERNAL_EXAHYPE_PATH, OPT_KERNEL_PATH_PREFIX));
  }
  
  //Remove all the files and subdir inside the given directory
  private static void cleanDirectory(final Path directory) throws IOException
  {
    if (Files.exists(directory))
    {
      Files.walkFileTree(directory, new SimpleFileVisitor<Path>()
      {
        @Override
        public FileVisitResult visitFile(Path path, BasicFileAttributes basicFileAttributes) throws IOException
        {
          Files.delete(path);
          return FileVisitResult.CONTINUE;
        }
        
        @Override
        public FileVisitResult postVisitDirectory(Path loc_dir, IOException ioException) throws IOException
        {
          if(!loc_dir.equals(directory))
            Files.delete(loc_dir);
          
          return FileVisitResult.CONTINUE;
        }
      });
    }
  }

}
