package eu.exahype;

import eu.exahype.analysis.DepthFirstAdapter;
import eu.exahype.node.AProject;

public class DirectoryAndPathChecker extends DepthFirstAdapter {
  public Boolean valid = true;
  
  private java.io.File peanoPath;
  private java.io.File tarchPath;
  private java.io.File multiscalelinkedcellPath;
  private java.io.File exahypePath;
  private java.io.File outputDirectory;
  
  public void inAProject(AProject node) {
	peanoPath                = new java.io.File(node.getPeanoPath().getText());
	tarchPath                = new java.io.File(node.getTarchPath().getText());
	multiscalelinkedcellPath = new java.io.File(node.getMultiscalelinkedcellPath().getText());
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
