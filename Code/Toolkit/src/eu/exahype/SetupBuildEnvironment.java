package eu.exahype;

import eu.exahype.analysis.DepthFirstAdapter;
import eu.exahype.node.AProject;
import eu.exahype.node.ASharedMemory;
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
      System.err.println( "ERROR: please set dimension to 2 instead of " + node.getDimension() + " if you specify a three-dimensional offset");
      valid = false;
	}
	else {
      System.out.print ("2d experiment ... ok\n" );      
	}

    try {
      _writer.write("PROJECT_CFLAGS+=-DDim2\n" );
    }
  	catch (Exception exc) {
      System.err.println( "ERROR: " + exc.toString() );
      valid = false;
  	}
  }

  @Override
  public void inAThreeDimensionalComputationalDomain(AThreeDimensionalComputationalDomain node) {
    if (!node.getDimension().toString().trim().equals("3")) {
	  System.err.println( "ERROR: please set dimension to 3 instead of " + node.getDimension() + " if you specify a three-dimensional offset");
	  valid = false;
	}
	else {
      System.out.print ("3d experiment ... ok\n" );      
    }

	try {
	  _writer.write("PROJECT_CFLAGS+=-DDim3\n" );
    }
  	catch (Exception exc) {
      System.err.println( "ERROR: " + exc.toString() );
      valid = false;
  	}
  }
  
  @Override
  public void inASharedMemory(ASharedMemory node) {
	try {
	  _writer.write("SHAREDMEM=TBB\n" );
      System.out.print ("shared memory ... TBB (switch to OpenMP manually as indicated below)\n" );   
      
      if (!System.getenv().containsKey("TBB_INC")) {
        System.out.print ("WARNING: environment variable TBB_INC not set but required if code is built with TBB\n" );      
      }
      if (!System.getenv().containsKey("TBB_SHLIB")) {
          System.out.print ("WARNING: environment variable TBB_SHLIB not set but required if code is built with TBB\n" );      
        }
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
      _writer.write("SHAREDMEMORYORACLES_PATH=" + _directoryAndPathChecker.sharedMemoryOraclesPath.getAbsolutePath() + "\n");
      _writer.write("EXAHYPE_PATH=" + _directoryAndPathChecker.exahypePath.getAbsolutePath() + "\n");
      _writer.write("PROJECT_PATH=" + _directoryAndPathChecker.outputDirectory.getAbsolutePath()  + "\n" );
      _writer.write("EXECUTABLE=ExaHyPE-" + node.getName() + "\n");
      _writer.write("\n\n");
      
      String architecture = node.getArchitecture().toString().trim().toLowerCase();
      
      switch(architecture) {
          case "wsm": 
              _writer.write("PROJECT_CFLAGS+=-DALIGNMENT=16");
              break;
          case "snb":
              _writer.write("PROJECT_CFLAGS+=-DALIGNMENT=32");
              break; 
          case "hsw": 
              _writer.write("PROJECT_CFLAGS+=-DALIGNMENT=32");
              break;
          case "knc":
              _writer.write("PROJECT_CFLAGS+=-DALIGNMENT=64");
              break;
          case "knl":
              _writer.write("PROJECT_CFLAGS+=-DALIGNMENT=64");
              break;
          default: 
              // noarch or unknown
              _writer.write("PROJECT_CFLAGS+=-DALIGNMENT=16");
              break;
      }
      _writer.write("\n");
      
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
      System.out.print ("  export CC=gcc  \t\t\tSelect GNU compiler\n");
      System.out.print ("  export CC=Intel\t\t\tSelect Intel compiler (default)\n");
      System.out.print ("\n");
      System.out.print ("  export MODE=Debug\t\t\tBuild debug version of code\n");
      System.out.print ("  export MODE=Asserts\t\t\tBuild release version of code that is augmented with assertions\n");
      System.out.print ("  export MODE=Profile\t\t\tBuild release version of code that produces profiling information\n");
      System.out.print ("  export MODE=Release\t\t\tBuild release version of code (default)\n");
      System.out.print ("\n");
      System.out.print ("  export SHAREDMEM=TBB\t\t\tUse Intel's Threading Building Blocks (TBB) for shared memory parallelisation\n");
      System.out.print ("  export SHAREDMEM=OMP\t\t\tUse OpenMP for shared memory parallelisation\n");
      System.out.print ("  export SHAREDMEM=\t\t\tDo not use shared memory (default if not indicated otherwise by \"shared memory ...\" message above)\n");
      System.out.print ("\n");
      System.out.print ("  export TBB_INC=-I...\t\t\tIndicate where to find TBB headers (only required if SHAREDMEM=TBB). Please add -I (Linux) prefix to path\n");
      System.out.print ("  export TBB_SHLIB=\"-L... -ltbb\"\tIndicate where to find TBB's shared libraries (only required if SHAREDMEM=TBB). Variable has to comprise both search path and library name\n");

      _writer.close();
	} 
	catch (Exception exc) {
      System.err.println( "ERROR: " + exc.toString() );
      valid = false;
	}
  }
}
