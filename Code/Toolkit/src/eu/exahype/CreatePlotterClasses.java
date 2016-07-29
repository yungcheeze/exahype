package eu.exahype;

import java.io.IOException;

import eu.exahype.analysis.DepthFirstAdapter;
import eu.exahype.node.AAderdgSolver;
import eu.exahype.node.AFiniteVolumesSolver;
import eu.exahype.node.AProfiling;
import eu.exahype.node.AProject;
import eu.exahype.node.APlotSolution;
import eu.exahype.node.AComputationalDomain;


public class CreatePlotterClasses extends DepthFirstAdapter {
  public Boolean valid = true;

  private DirectoryAndPathChecker _directoryAndPathChecker;

  private String                  _projectName;

  private int                     _dimensions;
  
  private String                  _solverName;
  
  private int                     _plotterCounter;

  private boolean                 _enableProfiler;

  public CreatePlotterClasses(DirectoryAndPathChecker directoryAndPathChecker) {
    _directoryAndPathChecker = directoryAndPathChecker;
  }

  @Override
  public void inAProject(AProject node) {
    _projectName     = node.getName().toString().trim();
  }

  @Override
  public void inAAderdgSolver(AAderdgSolver node) {
    _solverName     = node.getName().toString().trim();
    _plotterCounter = -1;
  }

  @Override
  public void inAFiniteVolumesSolver(AFiniteVolumesSolver node) {
    _solverName     = node.getName().toString().trim();
    _plotterCounter = -1;
  }

  @Override
  public void inAPlotSolution(APlotSolution node) {
    _plotterCounter++;
    try {
      java.io.File header         = new java.io.File(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/" + _solverName + "_Plotter" + Integer.toString(_plotterCounter) + ".h");
      java.io.File implementation = new java.io.File(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/" + _solverName + "_Plotter" + Integer.toString(_plotterCounter) + ".cpp");

      java.io.BufferedWriter headerWriter         = null;
      java.io.BufferedWriter implementationWriter = null;
          
      if (header.exists()) {
        System.err.println( "WARNING: file " + header.getName() + " already exists. Is not overwritten" );	
      }
      else {
        headerWriter = new java.io.BufferedWriter(new java.io.FileWriter(header.getAbsoluteFile()));
      }

      if (implementation.exists()) {
        System.err.println( "WARNING: file " + implementation.getName() + " already exists. Is not overwritten" );	
      }
      else {
        implementationWriter = new java.io.BufferedWriter(new java.io.FileWriter(implementation.getAbsoluteFile()));
      }

/*          if (
            node.getIdentifier().getText().trim().equals( "cellwise" )
            &&
            headerWriter!=null
          ) {
            writeCellWiseHeader(headerWriter,node);
          }
*/
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      valid = false;
    }
  }
}
