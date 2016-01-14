package eu.exahype;

import eu.exahype.analysis.DepthFirstAdapter;
import eu.exahype.node.AProject;

public class SetupBuildEnvironment extends DepthFirstAdapter {
  public Boolean valid = true;
  
  private DirectoryAndPathChecker   _directoryAndPathChecker;

  
  public SetupBuildEnvironment(DirectoryAndPathChecker  directoryAndPathChecker) {
	_directoryAndPathChecker = directoryAndPathChecker;
  }
  
  public void inAProject(AProject node) {
	if (!node.getDimensions().toString().trim().equals("2") && !node.getDimensions().toString().trim().equals("3")  ) {
      System.err.println( "ERROR: please set dimension to 2 or 3 instead of " + node.getDimensions() );
      valid = false;
	}
	try {
      java.io.File logFile = new java.io.File(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/Makefile");

      // @todo DIM has to be set here
      
      java.io.BufferedWriter writer = new java.io.BufferedWriter(new java.io.FileWriter(logFile));
      writer.write("PEANO_PATH=" + _directoryAndPathChecker.peanoPath.getAbsolutePath() + "\n");
      writer.write("TARCH_PATH=" + _directoryAndPathChecker.tarchPath.getAbsolutePath() + "\n");
      writer.write("MULTISCALELINKEDCELL_PATH=" + _directoryAndPathChecker.multiscalelinkedcellPath.getAbsolutePath() + "\n");
      writer.write("EXAHYPE_PATH=" + _directoryAndPathChecker.exahypePath.getAbsolutePath() + "\n");
      writer.write("PROJECT_PATH=" + _directoryAndPathChecker.outputDirectory.getAbsolutePath()  + "\n" );
      writer.write("DIM=-DDim" + node.getDimensions().toString().trim()  + "\n" );
      writer.write("\n\n\n\n");
      writer.write("-include " + _directoryAndPathChecker.exahypePath.getAbsolutePath() + "/Makefile\n");
      writer.write("\n\n\n\n");
      writer.write("all: \n");
      writer.write("\t@echo " + node.getName() + "\n");
      writer.write("\t@echo =================\n");
      writer.write("\t@echo An ExaHyPE solver\n");

      System.out.print ("store pathes in makefile ... ok" );      
	  System.out.println("\n\n\n\n");
      System.out.print ("please change into directory " + _directoryAndPathChecker.peanoPath.getAbsolutePath() + " and type make \n");
      System.out.print ("ensure that you set all environment variables before:\n");
      System.out.print ("  export CC=Intel\t\tSelect Intel compiler (default)\n");
      System.out.print ("  export CC=gcc  \t\tSelect GNU compiler\n");
      System.out.print ("\n");
      System.out.print ("  export MODE=Debug\t\tBuild debug version of code(default)\n");
      System.out.print ("  export MODE=Profile\t\tBuild release version of code that produces profiling information\n");
      System.out.print ("  export MODE=Release\t\tBuild release version of code\n");
      
      writer.close();
	} 
	catch (Exception exc) {
      System.err.println( "ERROR: " + exc.toString() );
      valid = false;
	}
  }
}
