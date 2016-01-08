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
      writer.write("Hello world!");
      
      writer.close();
	} 
	catch (Exception exc) {
      System.err.println( "ERROR: " + exc.toString() );
      valid = false;
	}
  }
}
