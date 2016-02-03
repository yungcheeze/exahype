package eu.exahype;

import eu.exahype.analysis.DepthFirstAdapter;
import eu.exahype.node.AProject;
import eu.exahype.node.APaths;

public class DirectoryAndPathChecker extends DepthFirstAdapter {
  public Boolean valid = true;
  
  protected java.io.File peanoPath;
  protected java.io.File tarchPath;
  protected java.io.File multiscalelinkedcellPath;
  protected java.io.File exahypePath;
  protected java.io.File outputDirectory;
  protected java.io.File sharedMemoryOraclesPath;

  @Override
  public void inAPaths(APaths node) {
	peanoPath                = new java.io.File(node.getPeanoPath().getText());
	tarchPath                = new java.io.File(node.getTarchPath().getText());
	multiscalelinkedcellPath = new java.io.File(node.getMultiscalelinkedcellPath().getText());
	sharedMemoryOraclesPath  = new java.io.File(node.getSharedmemoryoraclesPath().getText());
	exahypePath              = new java.io.File(node.getExahypePath().getText());
	outputDirectory          = new java.io.File(node.getOutputDirectory().getText());

	System.out.print ("Peano kernel path: " + peanoPath.getAbsolutePath() );
	if ( peanoPath.isDirectory() ) {
	  System.out.println( " ... ok" );
	}
	else {
	  System.out.println( " ... not found" );
	  valid = false;
	}
		
	System.out.print ("Peano tarch path: " + tarchPath.getAbsolutePath() );
	if ( tarchPath.isDirectory() ) {
	  System.out.println( " ... ok" );
	}
	else {
	  System.out.println( " ... not found" );
	  valid = false;
	}
		
	System.out.print ("multiscalelinkecell (Peano toolbox) path: " + multiscalelinkedcellPath.getAbsolutePath() );
	if ( multiscalelinkedcellPath.isDirectory() ) {
	  System.out.println( " ... ok" );
	}
	else {
	  System.out.println( " ... not found" );
	  valid = false;
	}
	
    System.out.print ("sharedmemoryoracles (Peano toolbox) path: " + sharedMemoryOraclesPath.getAbsolutePath() );
    if ( multiscalelinkedcellPath.isDirectory() ) {
      System.out.println( " ... ok" );
    }
    else {
      System.out.println( " ... not found" );
      valid = false;
    }
	
  };
  
  
  @Override
  public void outAProject(AProject node) {
	System.out.print ("ExaHyPE path: " + exahypePath.getAbsolutePath() );
	if ( exahypePath.isDirectory() ) {
      System.out.println( " ... ok" );
	}
	else {
	  System.out.println( " ... not found" );
	  valid = false;
	}
		
	System.out.print ("output directory: " + outputDirectory.getAbsolutePath() );
	if ( outputDirectory.isDirectory() ) {
      System.out.println( " ... does exist (will not be overwritten)" );
	}
	else {
      boolean createdDirectories = outputDirectory.mkdirs();
	  
      if (createdDirectories) {
        System.out.println( " ... created" );
      }
      else {
        System.out.println( " ... not found and could not be created" );
    	valid = false;
      }
	}
  }
}
