package eu.exahype.solvers;

import java.io.BufferedWriter;
import java.io.IOException;

public class GenericFluxesADER_DG implements Solver {
  public static final String Identifier = "generic::fluxes";
  protected boolean _hasConstants;
  protected boolean _isLinear;
  protected boolean _isFortran;
  
  private static final String fortranReference = "//************************************************* \n//for FORTRAN kernels the fluxes and eigenvalues \n//have to be implemented in the file ./PDE.f90 and ./typesDef.f90 \n//\n//You have to create these files yourself\n//and follow the sample applications in the wiki\n//************************************************* \n";

  public GenericFluxesADER_DG(int dimensions, int numberOfUnknowns, int numberOfParameters,
      int order, boolean enableProfiler, boolean hasConstants, boolean isLinear, boolean isFortran) {
    _dimensions = dimensions;
    _numberOfUnknowns = numberOfUnknowns;
    _numberOfParameters = numberOfParameters;
    _order = order;
    _enableProfiler = enableProfiler;
    _hasConstants   = hasConstants;
    _isLinear   = isLinear;
    _isFortran   = isFortran;
    
  }

  @Override
  public final void writeHeader(BufferedWriter writer, String solverName, String projectName)
      throws IOException {
     IncludeOnceHelper ifndef = new IncludeOnceHelper(writer, solverName+"_CLASS_HEADER");
     ifndef.open();
     Helpers.writeMinimalADERDGSolverHeader(solverName, writer, projectName, _hasConstants, _order, _dimensions, _numberOfUnknowns, _numberOfParameters);

    writer.write("\n");
    writer.write("    void init(std::vector<std::string>& cmdlineargs"+(_hasConstants ? ", exahype::Parser::ParserView& constants" : "")+");\n");
    writer.write("    void eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda);\n");
    writer.write("    void flux(const double* const Q, double** F);\n");
    writer.write("    void source(const double* const Q, double* S);\n");
    writer.write("    void boundaryValues(const double* const x,const double t, const double dt, const int faceIndex, const int normalNonZero, const double * const fluxIn, const double* const stateIn, double *fluxOut, double* stateOut);\n");
    writer.write(
        "    void adjustedSolutionValues(const double* const x,const double w,const double t,const double dt,double* Q);\n");
    writer.write(
        "    void ncp(const double* const Q, const double* const gradQ, double* BgradQ);\n");
    writer.write(
        "    void matrixb(const double* const Q, const int normalNonZero, double* Bn);\n");

    writer.write("};\n\n\n");
    ifndef.close();
  }

  @Override
  public final void writeGeneratedImplementation(BufferedWriter writer, String solverName,
      String projectName) throws IOException {

    writer.write("// ==============================================\n");
    writer.write("// Please do not change the implementations below\n");
    writer.write("// =============================---==============\n");
    writer.write("#include \"" + solverName + ".h\"\n");
    writer.write("#include \"kernels/aderdg/generic/Kernels.h\"\n");
    writer.write("\n\n");

    // constructor
    if (_hasConstants) {
      writer.write(projectName + "::" + solverName + "::" + solverName + "(double maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping, std::unique_ptr<exahype::profilers::Profiler> profiler, std::vector<std::string>& cmdlineargs, exahype::Parser::ParserView constants):\n");
    }
    else {
      writer.write(projectName + "::" + solverName + "::" + solverName + "(double maximumMeshSize, exahype::solvers::Solver::TimeStepping timeStepping, std::unique_ptr<exahype::profilers::Profiler> profiler, std::vector<std::string>& cmdlineargs):\n");
    }


    writer.write("  exahype::solvers::ADERDGSolver("
        + "\""+solverName+"\", nVar /* numberOfUnknowns */, "
        + "nParams /* numberOfParameters */, order + 1 "
        + " /* nodesPerCoordinateAxis */, maximumMeshSize, timeStepping, " +
        "std::move(profiler)) {\n");
    if(_hasConstants) {
       writer.write("  init(cmdlineargs, constants);\n");
    } else {
       writer.write("  init(cmdlineargs);\n");
       writer.write("  // PS: If you miss access to user constants here, enable them in the toolkit\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");

    String solverType = "<" + solverName + ">";
    String languageNamespace = _isFortran ? "fortran" : "c";

    writer.write("void " + projectName + "::" + solverName
        + "::spaceTimePredictor(double* lQhbnd,double* lFhbnd,double** tempSpaceTimeUnknowns,double** tempSpaceTimeFluxUnknowns,double*  tempUnknowns,double*  tempFluxUnknowns,const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& dx,const double dt) {\n");
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"spaceTimePredictor\");\n");
    }
    writer.write("  kernels::aderdg::generic::" + languageNamespace
        + "::spaceTimePredictor" + (_isLinear ? "Linear" : "Nonlinear") + solverType
        + "( *this, lQhbnd, lFhbnd, tempSpaceTimeUnknowns, tempSpaceTimeFluxUnknowns, tempUnknowns, tempFluxUnknowns, luh, dx, dt " + (_isFortran ? ", nVar, nParams, order + 1" : "")
      + ");\n");
    if (_enableProfiler) {
      writer.write("  _profiler->stop(\"spaceTimePredictor\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::solutionUpdate(double* luh, const double* const lduh, const double dt) {\n");
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"solutionUpdate\");\n");
    }
    writer.write("  kernels::aderdg::generic::" + languageNamespace
        + "::solutionUpdate( luh, lduh, dt, getNumberOfVariables(), getNumberOfParameters(), getNodesPerCoordinateAxis() );\n");
    if (_enableProfiler) {
      writer.write("  _profiler->stop(\"solutionUpdate\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::volumeIntegral(double* lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"volumeIntegral\");\n");
    }
    writer.write("  kernels::aderdg::generic::" + languageNamespace
        + "::volumeIntegral" + (_isLinear ? "Linear" : "Nonlinear")
        + "( lduh, lFhi, dx, getNumberOfVariables(), getNumberOfParameters(), getNodesPerCoordinateAxis() );\n");
    if (_enableProfiler) {
      writer.write("  _profiler->stop(\"volumeIntegral\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::surfaceIntegral(double* lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"surfaceIntegral\");\n");
    }
    writer.write("  kernels::aderdg::generic::" + languageNamespace
        + "::surfaceIntegral" + (_isLinear ? "Linear" : "Nonlinear")
        + "( lduh, lFhbnd, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    if (_enableProfiler) {
      writer.write("  _profiler->stop(\"surfaceIntegral\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::riemannSolver(double* FL, double* FR, const double* const QL, const double* const QR, double* tempFaceUnknownsArray, double** tempStateSizedVectors, double** tempStateSizedSquareMatrices, const double dt, const int normalNonZeroIndex) {\n");
    
    writer.write("  assertion2(normalNonZeroIndex>=0,dt,normalNonZeroIndex);\n");
    writer.write("  assertion2(normalNonZeroIndex<DIMENSIONS,dt,normalNonZeroIndex);\n");
    
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"riemannSolver\");\n");
    }
    writer.write("  kernels::aderdg::generic::" + languageNamespace
            + "::riemannSolver" + (_isLinear ? "Linear" : "Nonlinear") + solverType
            + "( *this, FL, FR, QL, QR, tempFaceUnknownsArray, tempStateSizedVectors, tempStateSizedSquareMatrices, dt, normalNonZeroIndex" + (_isFortran ? ", nVar, nParams, order + 1" : "")
    + " );\n");
    if (_enableProfiler) {
      writer.write("  _profiler->stop(\"riemannSolver\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");

    // boundaryConditions
    writer.write("void " + projectName + "::" + solverName
            + "::boundaryConditions(double* fluxOut,double* stateOut,const double* const fluxIn,const double* const stateIn,const tarch::la::Vector<DIMENSIONS, double>& cellCentre,const tarch::la::Vector<DIMENSIONS,double>& cellSize,const double t,const double dt,const int faceIndex,const int normalNonZero) {\n");
    if (_enableProfiler) {
        writer.write("  _profiler->start(\"boundaryConditions\");\n");
    }
    //ToDo only available as c++ implementation, reference it in Fortran namespace
    //writer.write("  kernels::aderdg::generic::" + languageNamespace
    writer.write("  kernels::aderdg::generic::c" 
            + "::boundaryConditions" + solverType
            + "( *this, fluxOut, stateOut, fluxIn, stateIn, cellCentre, cellSize, t, dt, faceIndex, normalNonZero );\n");
    if (_enableProfiler) {
        writer.write("  _profiler->stop(\"boundaryConditions\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
    
    writer.write("double " + projectName + "::" + solverName
        + "::stableTimeStepSize(const double* const luh,double* tempEigenvalues,const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"stableTimeStepSize\");\n");
    }
    writer.write("  double d = kernels::aderdg::generic::" + languageNamespace
        + "::stableTimeStepSize" + solverType
        + "( *this, luh, tempEigenvalues, dx" + (_isFortran ? ", nVar, order + 1" : "")
        + " );\n");
    if (_enableProfiler) {
      writer.write("  _profiler->stop(\"stableTimeStepSize\");\n");
    }
    writer.write("  return d;\n");
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::solutionAdjustment(double *luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {\n");
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"solutionAdjustment\");\n");
    }
    writer.write("  kernels::aderdg::generic::" + languageNamespace
        + "::solutionAdjustment" + solverType
        + "( *this, luh, center, dx, t, dt " + (_isFortran ? ", nVar, order + 1" : "")
        + ");\n");
    if (_enableProfiler) {
      writer.write("  _profiler->stop(\"solutionAdjustment\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::faceUnknownsProlongation(double* lQhbndFine,double* lFhbndFine,const double* lQhbndCoarse,const double* lFhbndCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) {\n");
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"faceUnknownsProlongation\");\n");
    }
    //ToDo only available as c++ implementation, reference it in Fortran namespace
    //writer.write("  kernels::aderdg::generic::" + languageNamespace
    writer.write("  kernels::aderdg::generic::c"
        + "::faceUnknownsProlongation( lQhbndFine, lFhbndFine, lQhbndCoarse, lFhbndCoarse, coarseGridLevel, fineGridLevel, subfaceIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    if (_enableProfiler) {
      writer.write("  _profiler->stop(\"faceUnknownsProlongation\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::faceUnknownsRestriction(double* lQhbndCoarse,double* lFhbndCoarse,const double* lQhbndFine,const double* lFhbndFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex) {\n");
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"faceUnknownsRestriction\");\n");
    }
    //ToDo only available as c++ implementation, reference it in Fortran namespace
    //writer.write("  kernels::aderdg::generic::" + languageNamespace
    writer.write("  kernels::aderdg::generic::c" 
        + "::faceUnknownsRestriction( lQhbndCoarse, lFhbndCoarse, lQhbndFine, lFhbndFine, coarseGridLevel, fineGridLevel, subfaceIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    if (_enableProfiler) {
      writer.write("  _profiler->stop(\"faceUnknownsRestriction\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::volumeUnknownsProlongation(double* luhFine, const double* luhCoarse, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {\n");
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"volumeUnknownsProlongation\");\n");
    }
    //ToDo only available as c++ implementation, reference it in Fortran namespace
    //writer.write("  kernels::aderdg::generic::" + languageNamespace
    writer.write("  kernels::aderdg::generic::c"
        + "::volumeUnknownsProlongation( luhFine, luhCoarse, coarseGridLevel, fineGridLevel, subcellIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    if (_enableProfiler) {
      writer.write("  _profiler->stop(\"volumeUnknownsProlongation\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");

    writer.write("void " + projectName + "::" + solverName
        + "::volumeUnknownsRestriction(double* luhCoarse, const double* luhFine, const int coarseGridLevel, const int fineGridLevel, const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) {\n");
    if (_enableProfiler) {
      writer.write("  _profiler->start(\"volumeUnknownsRestriction\");\n");
    }
    //ToDo only available as c++ implementation, reference it in Fortran namespace
    //writer.write("  kernels::aderdg::generic::" + languageNamespace
    writer.write("  kernels::aderdg::generic::c"
        + "::volumeUnknownsRestriction( luhCoarse, luhFine, coarseGridLevel, fineGridLevel, subcellIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
    if (_enableProfiler) {
      writer.write("  _profiler->stop(\"volumeUnknownsRestriction\");\n");
    }
    writer.write("}\n");
    writer.write("\n\n\n");
  }

  @Override
  public final void writeUserImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException {
    Helpers.writeMinimalADERDGSolverUserImplementation(solverName, writer, projectName,
        _numberOfUnknowns, _numberOfParameters, _order, _hasConstants);

    int digits = String.valueOf(_numberOfUnknowns + _numberOfParameters).length();

    
    if (_isFortran){
      writer.write(fortranReference);
    }else{

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
    }
    writer.write("\n\n\n");

    if (!_isFortran){
      // source
      writer.write("void " + projectName + "::" + solverName + "::source(const double* const Q, double* S) {\n");
      writer.write("  // Number of variables = " + _numberOfUnknowns + " + " +  _numberOfParameters + "\n");
      writer.write("  // @todo Please implement\n");
      for (int i = 0; i < _numberOfUnknowns + _numberOfParameters; i++) {
        writer.write("  S[" + i + "] = 0.0;\n");
      }
      writer.write("}\n");
    }
    writer.write("\n\n\n");
      
    // boundary conditions
    writer.write("void " + projectName + "::" + solverName
            + "::boundaryValues(const double* const x,const double t, const double dt, const int faceIndex, const int normalNonZero, const double * const fluxIn, const double* const stateIn, double *fluxOut, double* stateOut) {\n");
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
    
    if (!_isFortran){
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
    }
    writer.write("\n\n\n");
    
    //initial conditions
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

      // ncp
    if (!_isFortran){
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
    }
    writer.write("\n\n\n");

    if (!_isFortran){
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
    }
    writer.write("\n\n\n");
  }

  protected int _dimensions;
  protected int _numberOfUnknowns;
  protected int _numberOfParameters;
  protected int _order;
  protected boolean _enableProfiler;
}
