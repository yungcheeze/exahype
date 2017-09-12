package eu.exahype.solvers;

public interface Solver {
 
  /**
   * @return true if the solver supports generation of Variables classes.
   */
  public boolean supportsVariables();
  
  public String getSolverName();
  
  public void writeHeader(java.io.BufferedWriter writer) throws java.io.IOException;
  
  public void writeUserImplementation(java.io.BufferedWriter writer) throws java.io.IOException;
  
  public void writeAbstractHeader(java.io.BufferedWriter writer) throws java.io.IOException;
  
  public void writeAbstractImplementation(java.io.BufferedWriter writer) throws java.io.IOException;
}
