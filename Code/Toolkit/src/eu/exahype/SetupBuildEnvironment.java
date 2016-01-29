package eu.exahype;

import eu.exahype.analysis.DepthFirstAdapter;
import eu.exahype.node.AProject;
import eu.exahype.node.ATwoDimensionalComputationalDomain;
import eu.exahype.node.AThreeDimensionalComputationalDomain;

public class SetupBuildEnvironment extends DepthFirstAdapter {
  public Boolean valid = true;
  
  private java.io.BufferedWriter    _writer; 
  
  private DirectoryAndPathChecker   _directoryAndPathChecker;

  
  public SetupBuildEnvironment(DirectoryAndPathChecker  directoryAndPathChecker) {
	_directoryAndPathChecker = directoryAndPathChecker;
  }
  
  @Override
  public void inATwoDimensionalComputationalDomain(ATwoDimensionalComputationalDomain node) {
	if (!node.getDimension().toString().trim().equals("2")) {
      System.err.println( "WARNING: please set dimension to 2 instead of " + node.getDimension() + " if you specify a three-dimensional offset");
	}
	else {
      System.out.print ("2d experiment ... ok\n" );      
	}

    try {
      _writer.write("DIM=-DDim2\n" );
    }
  	catch (Exception exc) {
      System.err.println( "ERROR: " + exc.toString() );
      valid = false;
  	}
  }

  @Override
  public void inAThreeDimensionalComputationalDomain(AThreeDimensionalComputationalDomain node) {
    if (!node.getDimension().toString().trim().equals("3")) {
	  System.err.println( "WARNING: please set dimension to 3 instead of " + node.getDimension() + " if you specify a three-dimensional offset");
	}
	else {
      System.out.print ("3d experiment ... ok\n" );      
    }

	try {
	  _writer.write("DIM=-DDim3\n" );
    }
  	catch (Exception exc) {
      System.err.println( "ERROR: " + exc.toString() );
      valid = false;
  	}
  }
  
  @Override
  public void inAProject(AProject node) {
	try {
      java.io.File logFile = new java.io.File(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/Makefile");
      
      _writer = new java.io.BufferedWriter(new java.io.FileWriter(logFile));
      _writer.write("PEANO_PATH=" + _directoryAndPathChecker.peanoPath.getAbsolutePath() + "\n");
      _writer.write("TARCH_PATH=" + _directoryAndPathChecker.tarchPath.getAbsolutePath() + "\n");
      _writer.write("MULTISCALELINKEDCELL_PATH=" + _directoryAndPathChecker.multiscalelinkedcellPath.getAbsolutePath() + "\n");
      _writer.write("EXAHYPE_PATH=" + _directoryAndPathChecker.exahypePath.getAbsolutePath() + "\n");
      _writer.write("PROJECT_PATH=" + _directoryAndPathChecker.outputDirectory.getAbsolutePath()  + "\n" );
      _writer.write("EXECUTABLE=ExaHyPE-" + node.getName() + "\n");
      _writer.write("\n\n");
	} 
	catch (Exception exc) {
      System.err.println( "ERROR: " + exc.toString() );
      valid = false;
	}
  }
  
  @Override
  public void outAProject(AProject node) {
	try {
      _writer.write("\n\n");
      _writer.write("-include " + _directoryAndPathChecker.exahypePath.getAbsolutePath() + "/Makefile\n");
      _writer.write("\n\n\n\n");
      _writer.write("all: \n");
      _writer.write("\t@echo " + node.getName() + "\n");
      _writer.write("\t@echo =================\n");
      _writer.write("\t@echo An ExaHyPE solver\n");

      System.out.print ("store pathes and default settings in makefile ... ok" );      
	  System.out.println("\n\n\n\n");
      System.out.print ("please change into directory " + _directoryAndPathChecker.outputDirectory.getAbsolutePath() + " and type make \n");
      System.out.print ("ensure that you set all environment variables before:\n");
      System.out.print ("  export CC=gcc  \t\tSelect GNU compiler\n");
      System.out.print ("  export CC=Intel\t\tSelect Intel compiler (default)\n");
      System.out.print ("\n");
      System.out.print ("  export MODE=Debug\t\tBuild debug version of code\n");
      System.out.print ("  export MODE=Asserts\t\tBuild release version of code that is augmented with assertions\n");
      System.out.print ("  export MODE=Profile\t\tBuild release version of code that produces profiling information\n");
      System.out.print ("  export MODE=Release\t\tBuild release version of code (default)\n");
      
      _writer.close();
	} 
	catch (Exception exc) {
      System.err.println( "ERROR: " + exc.toString() );
      valid = false;
	}
  }
}
