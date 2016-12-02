package eu.exahype.solvers;

import eu.exahype.IOUtils;

public class LimitingADERDG implements Solver {
  public static final String Identifier = "generic::fluxes";

  private int _dimensions;
  private int _numberOfVariables;
  private int _numberOfParameters;
  private int _order;
//  private int _patchSize;
  private boolean _enableProfiler;
  private boolean _hasConstants;
  private boolean _isADERLinear;
  private boolean _isADERFortran;
  private boolean _isFVFortran;

  public LimitingADERDG(int dimensions, int numberOfVariables, int numberOfParameters,
      int order, boolean enableProfiler, boolean hasConstants, boolean isADERLinear, boolean isADERFortran, boolean isFVFortran) {
    _dimensions         = dimensions;
    _numberOfVariables  = numberOfVariables;
    _numberOfParameters = numberOfParameters;
    _order              = order;
//    _patchSize = patchSize;
    _enableProfiler     = enableProfiler;
    _hasConstants       = hasConstants;
    _isADERLinear       = isADERLinear;
    _isADERFortran      = isADERFortran;
    _isFVFortran        = isFVFortran;
  }

  public void writeHeader(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    
  }

  public void writeGeneratedImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException {
    
  }
  
  public void writeUserImplementation(java.io.BufferedWriter writer, String solverName,
      String projectName) throws java.io.IOException {
    
  }

  public void writeUserPDE(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    // @todo Implement
    System.err.println("Not implemented yet.\n");
  }


  public void writeTypesDef(java.io.BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    // @todo Implement
    System.err.println("Not implemented yet.\n");
  }
}
