package eu.exahype;

import eu.exahype.analysis.DepthFirstAdapter;
import eu.exahype.node.AProject;
import eu.exahype.node.PSolver;

public class GenerateSolverRegistration extends DepthFirstAdapter {
  public Boolean valid = true;

  private java.io.BufferedWriter _writer;
  private java.io.StringWriter _methodBodyWriter;

  private DirectoryAndPathChecker _directoryAndPathChecker;

  private int _kernelNumber;
  private int _plotterNumber;
  
  private String _solverName;

  private String _projectName;

  public GenerateSolverRegistration(DirectoryAndPathChecker directoryAndPathChecker) {
    _directoryAndPathChecker = directoryAndPathChecker;
    _kernelNumber = 0;
  }

  @Override
  public void inAProject(AProject node) {
    _projectName = node.getName().toString().trim();

    try {
      java.io.File logFile = new java.io.File(
          _directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/KernelCalls.cpp");

      _writer = new java.io.BufferedWriter(new java.io.FileWriter(logFile));
      _methodBodyWriter = new java.io.StringWriter();

      _writer.write("// This file is generated by the ExaHyPE toolkit.\n");
      _writer.write("// Please do not modify - it will be overwritten by the next\n");
      _writer.write("// ExaHyPE toolkit call.\n");
      _writer.write("// \n");
      _writer.write("// ========================\n");
      _writer.write("//   www.exahype.eu\n");
      _writer.write("// ========================\n\n");
      _writer.write("#include <sstream>\n\n");

      _writer.write("#include \"exahype/plotters/Plotter.h\"\n");
      _writer.write("#include \"exahype/profilers/ProfilerFactory.h\"\n");
      _writer.write("#include \"exahype/solvers/Solver.h\"\n");
      _writer.write("#include \"kernels/KernelCalls.h\"\n\n");

      _writer.write("#include \"kernels/GaussLegendreQuadrature.h\"\n");
      _writer.write("#include \"kernels/GaussLobattoQuadrature.h\"\n");
	  _writer.write("#include \"kernels/LimiterProjectionMatrices.h\"\n");
      _writer.write("#include \"kernels/DGMatrices.h\"\n");
      _writer.write("#include \"kernels/DGBasisFunctions.h\"\n\n");

      _methodBodyWriter.write("void kernels::initSolvers(exahype::Parser& parser) {\n");
      if (node.getSolver().size() == 0) {
        System.out.println("no solvers specified - create empty kernel calls ... ok");
      }
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      valid = false;
    }
  }
  
  private void writeProfilerCreation() {
      _methodBodyWriter.write("  std::string profiler_identifier = parser.getProfilerIdentifier();\n");
      _methodBodyWriter.write("  std::string metrics_identifier_list = parser.getMetricsIdentifierList();\n\n");

      _methodBodyWriter.write(
    		  "  assertion1(metrics_identifier_list.find_first_of(\"{\") == 0,\n"+
    		  "           metrics_identifier_list);\n");

      _methodBodyWriter.write(
    		  "  assertion1(metrics_identifier_list.find_last_of(\"}\") ==\n"+
    	      "                 metrics_identifier_list.size() - 1,\n"+
    		  "             metrics_identifier_list);\n\n");

      _methodBodyWriter.write("  // Split \"{metric1,metric2...}\" into {\"metric1\", \"metric2\", ...}\n");
      _methodBodyWriter.write("  std::vector<std::string> metrics_vector;\n");
      _methodBodyWriter.write("  std::stringstream ss;\n");
      _methodBodyWriter.write("  ss << metrics_identifier_list.substr(1, metrics_identifier_list.size() - 2);\n");
      _methodBodyWriter.write("  std::string metric;\n");
      _methodBodyWriter.write(
          "  while (std::getline(ss, metric, ',')) {\n"+
          "    metrics_vector.emplace_back(std::move(metric));\n"+
          "  }\n\n");

      _methodBodyWriter.write("  // Create profiler\n");
      _methodBodyWriter.write(
    		  "  auto profiler = exahype::profilers::ProfilerFactory::getInstance().create(\n"+
    		  "    profiler_identifier, metrics_vector);\n\n");
  }

  @Override
  public void inAAderdgSolver(eu.exahype.node.AAderdgSolver node) {
    try {
      _solverName = node.getName().toString().trim();
      int order         = Integer.parseInt(node.getOrder().getText());

      _writer.write("#include \"" + _solverName + ".h\"\n");

      _methodBodyWriter.write("  {\n");
      
      writeProfilerCreation();
      
      _methodBodyWriter.write("  // Create and register solver\n");
      _methodBodyWriter.write("  exahype::solvers::RegisteredSolvers.push_back( new " + _projectName +
    		                  "::" + _solverName + "("+order+"+1, parser.getMaximumMeshSize("+_kernelNumber+"), parser.getTimeStepping("+_kernelNumber+"), std::move(profiler)\n");
      if (node.getConstants()!=null) {
          _methodBodyWriter.write( "  , parser.getParserView(" +  _kernelNumber + ")\n");
        }
      _methodBodyWriter.write( "  ));\n");
      _methodBodyWriter.write("  parser.checkSolverConsistency("+_kernelNumber+");\n\n");
      _methodBodyWriter.write("  \n");
  
      _methodBodyWriter.write("  }\n");
      
      _kernelNumber++;
      _plotterNumber = 0;

      System.out.println("added creation of solver " + _solverName + " ... ok");
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      valid = false;
    }
  };

  @Override
  public void inAFiniteVolumesSolver(eu.exahype.node.AFiniteVolumesSolver node) {
    try {
      _solverName = node.getName().toString().trim();
      int patchSize     = Integer.parseInt(node.getPatchSize().getText());
      
      _writer.write("#include \"" + _solverName + ".h\"\n");

      _methodBodyWriter.write("  {\n"); // why do we need this?
      
      writeProfilerCreation();

      _methodBodyWriter.write("  // Create and register solver\n");
      _methodBodyWriter.write("  exahype::solvers::RegisteredSolvers.push_back( new " + _projectName +
    		                  "::" + _solverName + "("+patchSize+", parser.getMaximumMeshSize("+_kernelNumber+"), parser.getTimeStepping("+_kernelNumber+"), std::move(profiler)" );
      if (node.getConstants()!=null) {
        _methodBodyWriter.write( "  , parser.getParserView(" +  _kernelNumber + ")\n");
      }
      _methodBodyWriter.write( "  ));\n");
      _methodBodyWriter.write("  parser.checkSolverConsistency("+_kernelNumber+");\n\n");
      _methodBodyWriter.write("  \n");
      
      _methodBodyWriter.write("  }\n"); // why do we need this?
      
      _kernelNumber++;
      _plotterNumber = 0;

      System.out.println("added creation of solver " + _solverName + " ... ok");
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      valid = false;
    }
  };

  @Override
  public void inAPlotSolution(eu.exahype.node.APlotSolution node) {
    try {
      String plotterName = _solverName + "_Plotter" + Integer.toString(_plotterNumber);

      _writer.write("#include \"" + plotterName + ".h\"\n");

      _methodBodyWriter.write(
          "  exahype::plotters::RegisteredPlotters.push_back( new exahype::plotters::Plotter("
              + (_kernelNumber - 1) + "," + _plotterNumber + ",parser,new " + _projectName + "::" + plotterName + "(  *static_cast<" + _projectName + "::" + _solverName + "*>(exahype::solvers::RegisteredSolvers[" + (_kernelNumber-1) + "])) ));\n\n");
      _plotterNumber++;
      System.out.println("added plotter ... ok");
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      valid = false;
    }
  };

  @Override
  public void outAProject(AProject node) {
    try {
      _methodBodyWriter.write("\n");
      _methodBodyWriter.write(
          "  std::set<int> orders;\n"+
          "  for (const auto p : exahype::solvers::RegisteredSolvers) {\n"+
          "    orders.insert(p->getNodesPerCoordinateAxis()-1);\n"+
          "  }\n"+
          "  kernels::initGaussLegendreNodesAndWeights(orders);\n"+
          "  kernels::initGaussLobattoNodesAndWeights(orders);\n"+
          "  kernels::initLimiterProjectionMatrices(orders);\n"+
          "  kernels::initDGMatrices(orders);\n" +
          "  kernels::initBasisFunctions(orders);\n");
      _methodBodyWriter.write("}\n"); // close initSolvers(...)
      _methodBodyWriter.write("\n");
      _methodBodyWriter.write("\n");
      _methodBodyWriter.write("void kernels::finalise() {\n");
      _methodBodyWriter.write(
          "  std::set<int> orders;\n"+
          "  for (const auto p : exahype::solvers::RegisteredSolvers) {\n" +
          "    orders.insert(p->getNodesPerCoordinateAxis()-1);\n" +
          "  }\n" +
          "  kernels::freeGaussLegendreNodesAndWeights(orders);\n"+
          "  kernels::freeGaussLobattoNodesAndWeights(orders);\n"+
          "  kernels::freeLimiterProjectionMatrices(orders);\n"+
          "  kernels::freeDGMatrices(orders);\n"+
          "  kernels::freeBasisFunctions(orders);\n\n"+
          "  for (auto solver : exahype::solvers::RegisteredSolvers) {\n"+
          "    delete solver;\n"+
          "  }\n"+
          "  exahype::solvers::RegisteredSolvers.clear();\n\n"+
          "  for (auto plotter : exahype::plotters::RegisteredPlotters) {\n"+
          "    delete plotter;\n"+
          "  }\n"+
          "  exahype::plotters::RegisteredPlotters.clear();\n");
      _methodBodyWriter.write("}\n");
      _writer.write("\n");
      _writer.write("\n");
      _writer.write("\n");
      _writer.write(_methodBodyWriter.toString());
      _writer.write("\n");

      System.out.println("configured all solver solvers ... ok");

      _writer.write("\n\n");
      _writer.close();
    } catch (Exception exc) {
      System.err.println("ERROR: " + exc.toString());
      valid = false;
    }
  }
}
