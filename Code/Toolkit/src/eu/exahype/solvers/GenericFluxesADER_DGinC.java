package eu.exahype.solvers;

import java.io.BufferedWriter;
import java.io.IOException;

public abstract class GenericFluxesADER_DGinC extends GenericFluxesADER_DG {

  public GenericFluxesADER_DGinC(int dimensions, int numberOfUnknowns, int numberOfParameters,
      int order, boolean enableProfiler) {
    super(dimensions, numberOfUnknowns, numberOfParameters, order, enableProfiler);
  }

  @Override
  public final void writeUserImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException {
    Helpers.writeMinimalADERDGSolverUserImplementation(solverName, writer, projectName,
        _numberOfUnknowns, _numberOfParameters, _order);

    int digits = String.valueOf(_numberOfUnknowns + _numberOfParameters).length();

    // flux
      writer.write("void " + projectName + "::" + solverName
          + "::flux(const double* const Q, double** F) {\n");
    writer.write("  // Dimensions             = " + _dimensions + "\n");
    writer.write(
        "  // Number of variables    = " + Integer.toString(_numberOfUnknowns + _numberOfParameters)
            + " (#unknowns + #parameters)\n\n");
    writer.write("  double* f = F[0];\n");
    writer.write("  double* g = F[1];\n");
    if (_dimensions == 3) {
      writer.write("  double* h = F[2];\n");
    }
    writer.write("\n");
    writer.write("  // @todo Please implement\n");
    writer.write("  // f\n");
    for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
      writer.write("  f[" + String.format("%" + digits + "d", i) + "] = 0.0;\n");
    }
    writer.write("  // g\n");
    writer.write("  // @todo Please implement\n");
    for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
      writer.write("  g[" + String.format("%" + digits + "d", i) + "] = 0.0;\n");
    }
    if (_dimensions == 3) {
      writer.write("  // h\n");
      writer.write("  // @todo Please implement\n");
      for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
        writer.write("  h[" + String.format("%" + digits + "d", i) + "] = 0.0;\n");
      }
    }
    writer.write("}\n");
    writer.write("\n\n\n");

    // source
    writer.write("void " + projectName + "::" + solverName + "::source(const double* const Q, double* S) {\n");
    writer.write("  // Number of variables = " + _numberOfUnknowns + " + " +  _numberOfParameters + "\n");
    writer.write("  // @todo Please implement\n");
    for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
      writer.write("  S[" + i + "] = 0.0;\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
    
    // boundary conditions
    writer.write("void " + projectName + "::" + solverName
            + "::boundaryValues(const double* const x,const double t, const int faceIndex, const int normalNonZero, const double * const fluxIn, const double* const stateIn, double *fluxOut, double* stateOut) {\n");
    writer.write("  // Dimensions             = " + _dimensions + "\n");
    writer.write(
            "  // Number of variables    = " + Integer.toString(_numberOfUnknowns + _numberOfParameters)
            + " (#unknowns + #parameters)\n\n");
    writer.write("\n");
    writer.write("  // @todo Please implement\n");
    writer.write("  // fluxOut\n");
    for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
        writer.write("  fluxOut[" + String.format("%" + digits + "d", i) + "] = fluxIn[" + String.format("%" + digits + "d", i) + "];\n");
    }
    writer.write("  // stateOut\n");
    writer.write("  // @todo Please implement\n");
    for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
        writer.write("  stateOut[" + String.format("%" + digits + "d", i) + "] = stateIn[" + String.format("%" + digits + "d", i) + "];\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
    
    // eigenvalues
    writer.write("void " + projectName + "::" + solverName
        + "::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {\n");
    writer.write("  // Dimensions             = " + _dimensions + "\n");
    writer.write(
        "  // Number of variables    = " + Integer.toString(_numberOfUnknowns + _numberOfParameters)
            + " (#unknowns + #parameters)\n");
    writer.write("  // @todo Please implement\n");
    for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
      writer.write("  lambda[" + String.format("%" + digits + "d", i) + "] = 0.0;\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("bool " + projectName + "::" + solverName
        + "::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double> &center, const tarch::la::Vector<DIMENSIONS, double> &dx, double t) {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("  return false;\n");
    writer.write("}\n");
    writer.write("\n\n\n");
    writer.write("void " + projectName + "::" + solverName
        + "::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q) {\n");
    writer.write("  // Dimensions             = " + _dimensions + "\n");
    writer.write(
        "  // Number of variables    = " + Integer.toString(_numberOfUnknowns + _numberOfParameters)
            + " (#unknowns + #parameters)\n");
    writer.write("  // @todo Please implement\n");
    for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
      writer.write("  Q[" + String.format("%" + digits + "d", i) + "] = 0.0;\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");

    // refinement control
    writer.write("exahype::solvers::Solver::RefinementControl " + projectName + "::" + solverName
        + "::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, double t, const int level) {\n");
    writer.write("  // @todo Please implement\n");
    writer.write("  return exahype::solvers::Solver::RefinementControl::Keep;\n");
    writer.write("}\n");
    writer.write("\n\n\n");

    if (isLinear()) {
      // ncp
      writer.write("void " + projectName + "::" + solverName
          + "::ncp(const double* const Q, const double* const gradQ, double* BgradQ) {\n");
      writer.write("  // Dimensions             = " + _dimensions + "\n");
      writer.write("  // Number of variables    = "
          + Integer.toString(_numberOfUnknowns + _numberOfParameters)
          + " (#unknowns + #parameters)\n");
      writer.write("  // @todo Please implement\n");
      for (int i = 0; i < _dimensions * (_numberOfUnknowns + _numberOfParameters); i++) {
          writer.write("  BgradQ[" + i + "] = 0.0;\n");
      }
      writer.write("}\n");
      writer.write("\n\n\n");

      // matrixb
      writer.write("void " + projectName + "::" + solverName + "::matrixb(const double* const Q, const int normalNonZero, double* Bn) {\n");
      writer.write("  // Number of variables    = "
          + Integer.toString(_numberOfUnknowns + _numberOfParameters)
          + " (#unknowns + #parameters)\n");
      writer.write("  // @todo Please implement\n");
      for (int i = 0; i < (_numberOfUnknowns + _numberOfParameters) * (_numberOfUnknowns + _numberOfParameters); i++) {
        writer.write("Bn[" + i + "] = 0.0;\n");
      }
      writer.write("}\n");
      writer.write("\n\n\n");
    }
  }

  @Override
  public final void writeUserPDE(BufferedWriter writer, String solverName, String projectName) throws IOException {
    System.err.println("C-style kernels do not have a PDF.f90.\n");
  }

  @Override
  public final void writeTypesDef(BufferedWriter writer, String solverName, String projectName)
      throws IOException {
    System.err.println("C-style kernels do not have a typesDef.f90.\n");
  }

}
