package eu.exahype.solvers;

public interface Solver {
 
  /**
   * Configuration parameter: id of the options
   * (static final is implicit)
   */
  public String LINEAR_OPTION_ID      = "linear";
  public String NONLINEAR_OPTION_ID   = "nonlinear";
  public String FLUX_OPTION_ID        = "fluxes";
  public String SOURCE_OPTION_ID      = "sources";
  public String NCP_OPTION_ID         = "ncp";
  public String NO_TIME_AVG_OPTION_ID = "notimeavg";
 
  /**
   * @return true if the solver supports generation of Variables classes.
   */
  public boolean supportsVariables();
  
  public String getSolverName();
  
  public void writeHeader(java.io.BufferedWriter writer) throws java.io.IOException, IllegalArgumentException;
  
  public void writeUserImplementation(java.io.BufferedWriter writer) throws java.io.IOException, IllegalArgumentException;
  
  public void writeAbstractHeader(java.io.BufferedWriter writer) throws java.io.IOException, IllegalArgumentException;
  
  public void writeAbstractImplementation(java.io.BufferedWriter writer) throws java.io.IOException, IllegalArgumentException;

}
