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
	try {
      java.io.File logFile = new java.io.File(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/Makefile");

      java.io.BufferedWriter writer = new java.io.BufferedWriter(new java.io.FileWriter(logFile));
      writer.write("PEANO_PATH=" + _directoryAndPathChecker.peanoPath.getAbsolutePath() + "\n");
      writer.write("TARCH_PATH=" + _directoryAndPathChecker.tarchPath.getAbsolutePath() + "\n");
      writer.write("MULTISCALELINKEDCELL_PATH=" + _directoryAndPathChecker.multiscalelinkedcellPath.getAbsolutePath() + "\n");
      writer.write("EXAHYPE_PATH=" + _directoryAndPathChecker.exahypePath.getAbsolutePath() + "\n");
      writer.write("PROJECT_PATH=" + _directoryAndPathChecker.outputDirectory.getAbsolutePath() );
      writer.write("\n\n\n\n");
      writer.write("-include " + _directoryAndPathChecker.exahypePath.getAbsolutePath() + "/Makefile\n");
      writer.write("\n\n\n\n");
      writer.write("all: \n");
      writer.write("\t@echo " + node.getName() + "\n");
      writer.write("\t@echo =================\n");
      writer.write("\t@echo An ExaHyPE solver\n");

      System.out.print ("store pathes in makefile ... ok" );      
      
      writer.close();
	} 
	catch (Exception exc) {
      System.err.println( "ERROR: " + exc.toString() );
      valid = false;
	}
  }
}
