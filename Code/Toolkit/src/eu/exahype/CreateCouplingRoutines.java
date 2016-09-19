package eu.exahype;

import java.io.IOException;

import eu.exahype.analysis.DepthFirstAdapter;
import eu.exahype.node.AAderdgSolver;
import eu.exahype.node.ACoupleSolvers;
import eu.exahype.node.PSolver;
import eu.exahype.node.AFiniteVolumesSolver;
import eu.exahype.node.AProfiling;
import eu.exahype.node.AProject;
import eu.exahype.node.AComputationalDomain;


public class CreateCouplingRoutines extends DepthFirstAdapter {
  public Boolean valid = true;
  
  private DirectoryAndPathChecker _directoryAndPathChecker;
  
  private int       _numberOfSolvers;
  private String    _projectName;
  private AProject  _project;
  
  public CreateCouplingRoutines(DirectoryAndPathChecker directoryAndPathChecker) {
    _directoryAndPathChecker = directoryAndPathChecker;
  }


  @Override
  public void inAProject(AProject node) {
    _numberOfSolvers = node.getSolver().size();
    _projectName     = node.getName().getText().trim();
    _project         = node;
  }
  
  
  private void writeSolverIncludes(java.io.BufferedWriter writer) throws java.io.IOException {
    for (PSolver from: _project.getSolver()) {
      if (from instanceof eu.exahype.node.AAderdgSolver ) {
        eu.exahype.node.AAderdgSolver fromADERDG = (eu.exahype.node.AAderdgSolver)from;
        writer.write( "#include \"" + fromADERDG.getName().toString().trim() + ".h\"\n" );
      }
      if (from instanceof eu.exahype.node.AFiniteVolumesSolver ) {
        eu.exahype.node.AFiniteVolumesSolver fromFiniteVolumes = (eu.exahype.node.AFiniteVolumesSolver)from;
        writer.write( "#include \"" + fromFiniteVolumes.getName().toString().trim() + ".h\"\n" );
      }
    }
    writer.write( "\n\n\n" );
  }

  
  private void writeSolverForwardDeclarations(java.io.BufferedWriter writer) throws java.io.IOException {
    writer.write("namespace " + _projectName + "{\n");
    for (PSolver from: _project.getSolver()) {
      if (from instanceof eu.exahype.node.AAderdgSolver ) {
        eu.exahype.node.AAderdgSolver fromADERDG = (eu.exahype.node.AAderdgSolver)from;
        writer.write( "  class " + fromADERDG.getName().toString().trim() + ";\n" );
      }
      if (from instanceof eu.exahype.node.AFiniteVolumesSolver ) {
        eu.exahype.node.AFiniteVolumesSolver fromFiniteVolumes = (eu.exahype.node.AFiniteVolumesSolver)from;
        writer.write( "  class " + fromFiniteVolumes.getName().toString().trim() + ";\n" );
      }
    }
    writer.write("}\n\n\n");
  }

  
  private void writeCellWiseHeader(java.io.BufferedWriter headerWriter, ACoupleSolvers node) throws java.io.IOException {
    eu.exahype.solvers.Helpers.writeHeaderCopyright(headerWriter);
    headerWriter.write( "#include \"exahype/solvers/CellWiseCoupling.h\"\n" );
    headerWriter.write( "\n\n\n" );
    writeSolverForwardDeclarations(headerWriter);
    headerWriter.write("namespace " + _projectName + "{\n");
    headerWriter.write("  class " + node.getIdentifier().getText().trim() + ";\n");
    headerWriter.write("}\n\n\n");

    headerWriter.write("class " + _projectName + "::" + node.getIdentifier().getText().trim() + ": public exahype::solvers::CellWiseCoupling {\n");
    headerWriter.write("public:\n");
    headerWriter.write("  " + node.getIdentifier().getText().trim() + "(double time, double repeat); \n");
    headerWriter.write("\n\n");

    for (PSolver from: _project.getSolver()) 
    for (PSolver to:   _project.getSolver()) {
      headerWriter.write( "  /**\n" );
      headerWriter.write( "   * @todo Add your comment here\n" );
      headerWriter.write( "   *       This file is not overwritten if you rerun the ExaHyPE toolkit. If you \n" );
      headerWriter.write( "   *       add new solvers and thus need updated coupling routines, you have to \n" );
      headerWriter.write( "   *       delete this file manually before you rerun the toolkit. \n" );
      headerWriter.write( "   */\n" );
      headerWriter.write( "  void couple(\n" );
      if (from instanceof eu.exahype.node.AAderdgSolver ) {
        eu.exahype.node.AAderdgSolver fromADERDG = (eu.exahype.node.AAderdgSolver)from;
        headerWriter.write( "    " + fromADERDG.getName().toString().trim() + "&  fromSolver, \n" );
        //headerWriter.write( "    exahype.records.ADERDGCellDescription&  fromDescription, \n" );
        headerWriter.write( "    double*  fromluh, \n" );
        headerWriter.write( "    double*  fromlQhbnd, \n" );
        headerWriter.write( "    bool     fromHoldsValidMinMaxCondition, \n" );
      }
      if (from instanceof eu.exahype.node.AFiniteVolumesSolver ) {
        eu.exahype.node.AFiniteVolumesSolver fromFiniteVolumes = (eu.exahype.node.AFiniteVolumesSolver)from;
        headerWriter.write( "    " + fromFiniteVolumes.getName().toString().trim() + "&  fromSolver, \n" );
        //headerWriter.write( "    exahype.records.FiniteVolumesCellDescription&  fromDescription, \n" );
        headerWriter.write( "    double*  fromluh, \n" );
      }
      if (to instanceof eu.exahype.node.AAderdgSolver ) {
        eu.exahype.node.AAderdgSolver fromADERDG = (eu.exahype.node.AAderdgSolver)to;
        headerWriter.write( "    " + fromADERDG.getName().toString().trim() + "&  toSolver, \n" );
        //headerWriter.write( "    exahype.records.ADERDGCellDescription&  toDescription, \n" );
        headerWriter.write( "    double*  toluh, \n" );
        headerWriter.write( "    double*  tolQhbnd, \n" );
        headerWriter.write( "    bool     toHoldsValidMinMaxCondition\n" );
      }
      if (to instanceof eu.exahype.node.AFiniteVolumesSolver ) {
        eu.exahype.node.AFiniteVolumesSolver fromFiniteVolumes = (eu.exahype.node.AFiniteVolumesSolver)to;
        headerWriter.write( "    " + fromFiniteVolumes.getName().toString().trim() + "&  toSolver, \n" );
        //headerWriter.write( "    exahype.records.FiniteVolumesCellDescription&  toDescription, \n" );
        headerWriter.write( "    double*  toluh\n" );
      }
      headerWriter.write( "  );\n\n\n\n" );
    }
    headerWriter.write("};\n\n\n");
  }

  
  private void writeCellWiseImplementation(java.io.BufferedWriter implementationWriter, ACoupleSolvers node) throws java.io.IOException {
    implementationWriter.write( "#include \"" + node.getIdentifier().getText().trim() + ".h\"\n" );
    implementationWriter.write( "\n" );
    implementationWriter.write( "\n" );
    writeSolverIncludes(implementationWriter);
    implementationWriter.write( "\n" );

    implementationWriter.write( _projectName + "::" + node.getIdentifier().getText().trim() + "::" + node.getIdentifier().getText().trim() + "(double time, double repeat):\n" );
    implementationWriter.write( "  exahype::solvers::CellWiseCoupling(time,repeat) {\n" );
    implementationWriter.write( "  // @todo add your code here\n" );
    implementationWriter.write( "}" );
    implementationWriter.write( "\n" );
    implementationWriter.write( "\n" );
    implementationWriter.write( "\n" );
    
    AProject project = (AProject)(node.parent());
    for (PSolver from: project.getSolver()) 
    for (PSolver to:   project.getSolver()) {
      implementationWriter.write( "/**\n" );
      implementationWriter.write( " * This file is not overwritten if you rerun the ExaHyPE toolkit. If you \n" );
      implementationWriter.write( " * add new solvers and thus need updated coupling routines, you have to \n" );
      implementationWriter.write( " * delete this file manually before you rerun the toolkit. \n" );
      implementationWriter.write( " */\n" );
      implementationWriter.write( "void " + _projectName + "::" + node.getIdentifier().getText().trim() + "::couple(\n" );
      if (from instanceof eu.exahype.node.AAderdgSolver ) {
        eu.exahype.node.AAderdgSolver fromADERDG = (eu.exahype.node.AAderdgSolver)from;
        implementationWriter.write( "  " + fromADERDG.getName().toString().trim() + "&  fromSolver, \n" );
        //implementationWriter.write( "  exahype.records.ADERDGCellDescription&  fromDescription, \n" );
        implementationWriter.write( "  double*  fromluh, \n" );
        implementationWriter.write( "  double*  fromlQhbnd, \n" );
        implementationWriter.write( "  bool     fromHoldsValidMinMaxCondition, \n" );
      }
      if (from instanceof eu.exahype.node.AFiniteVolumesSolver ) {
        eu.exahype.node.AFiniteVolumesSolver fromFiniteVolumes = (eu.exahype.node.AFiniteVolumesSolver)from;
        implementationWriter.write( "  " + fromFiniteVolumes.getName().toString().trim() + "&  fromSolver, \n" );
        //implementationWriter.write( "  exahype.records.FiniteVolumesCellDescription&  fromDescription, \n" );
        implementationWriter.write( "  double*  fromluh, \n" );
      }
      if (to instanceof eu.exahype.node.AAderdgSolver ) {
        eu.exahype.node.AAderdgSolver fromADERDG = (eu.exahype.node.AAderdgSolver)to;
        implementationWriter.write( "  " + fromADERDG.getName().toString().trim() + "&  toSolver, \n" );
        //implementationWriter.write( "  exahype.records.ADERDGCellDescription&  toDescription, \n" );
        implementationWriter.write( "  double*  toluh, \n" );
        implementationWriter.write( "  double*  tolQhbnd, \n" );
        implementationWriter.write( "  bool     toHoldsValidMinMaxCondition\n" );
      }
      if (to instanceof eu.exahype.node.AFiniteVolumesSolver ) {
        eu.exahype.node.AFiniteVolumesSolver fromFiniteVolumes = (eu.exahype.node.AFiniteVolumesSolver)to;
        implementationWriter.write( "  " + fromFiniteVolumes.getName().toString().trim() + "&  toSolver, \n" );
        //implementationWriter.write( "  exahype.records.FiniteVolumesCellDescription&  toDescription, \n" );
        implementationWriter.write( "  double*  toluh\n" );
      }
      implementationWriter.write( "{\n" );
      implementationWriter.write( "  // @todo add your code here \n" );
      implementationWriter.write( "}\n\n\n\n" );
    }
  }

	    
  @Override
  public void inACoupleSolvers(ACoupleSolvers node) {
  	if (node.getType().getText().trim().equals( "cellwise" ) ) {
      System.out.println( "- generate cell-wise coupling layer" );
    }
  	else {
      System.err.println( "ERROR: coupling identifier " + node.getType().getText().trim() + " is not supported" );
      valid = false;
  	}

    if (_numberOfSolvers==1) {
      System.err.println( "WARNING: there is only one solver specified. Coupling might be irrelevant" );
    }

  	if (valid) {
      try {
        java.io.File header         = new java.io.File(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/" + node.getIdentifier().getText().trim() + ".h");
        java.io.File implementation = new java.io.File(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/" + node.getIdentifier().getText().trim() + ".cpp");

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

        if (
          node.getType().getText().trim().equals( "cellwise" )
          &&
          headerWriter!=null
        ) {
          writeCellWiseHeader(headerWriter,node);
        }

        if (
          node.getType().getText().trim().equals( "cellwise" )
          &&
          implementationWriter!=null
        ) {
          writeCellWiseImplementation(implementationWriter,node);
        }
        
        if (headerWriter!=null) {
          headerWriter.close();
        }
        if (implementationWriter!=null) {
          implementationWriter.close();
        }
      } catch (Exception exc) {
        System.err.println("ERROR: " + exc.toString());
        valid = false;
      }
  	}
  }
}
