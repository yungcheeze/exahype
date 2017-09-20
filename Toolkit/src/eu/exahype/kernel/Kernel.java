package eu.exahype.kernel;

import java.util.Set;
import java.util.stream.Collectors;

import eu.exahype.node.PIds;
import eu.exahype.node.AIds;
import eu.exahype.node.AIdentifierId;

import eu.exahype.node.AAderdgSolver;
import eu.exahype.node.ALimitingAderdgSolver;
import eu.exahype.node.PSolver;

public class Kernel {
  
   
  /**
   * Configuration parameter: id of the options
   */
  public static final String LINEAR_OPTION_ID      = "linear";
  public static final String NONLINEAR_OPTION_ID   = "nonlinear";
  public static final String FLUX_OPTION_ID        = "flux";
  public static final String SOURCE_OPTION_ID      = "source";
  public static final String NCP_OPTION_ID         = "ncp";
  public static final String NO_TIME_AVG_OPTION_ID = "notimeavg";
  public static final String POINTSOURCE_OPTION_ID = "pointsources";
  
  private Set<String> type;
  private Set<String> terms;
  private Set<String> optimization;
  
  public Kernel(PSolver solver) throws IllegalArgumentException {
    if(solver instanceof AAderdgSolver) {
      type = parseIds(((AAderdgSolver) solver).getKernelType());
      terms = parseIds(((AAderdgSolver) solver).getKernelTerms());
      optimization = parseIds(((AAderdgSolver) solver).getKernelOpt());
    } else if(solver instanceof ALimitingAderdgSolver) {
      type = parseIds(((ALimitingAderdgSolver) solver).getKernelType());
      terms = parseIds(((ALimitingAderdgSolver) solver).getKernelTerms());
      optimization = parseIds(((ALimitingAderdgSolver) solver).getKernelOpt());
    } else {
      throw new IllegalArgumentException("No kernel definition found");
    }
    
    validate();
  }
  
  //return null on error, use only after the program should already have failed with invalid kernel
  public static Kernel noExceptionContructor(PSolver solver) {
    try {
      return new Kernel(solver);
    } catch(IllegalArgumentException e) {
      return null;
    }
  }
  
  private static Set<String> parseIds(PIds idsRaw) {
    return ((AIds)idsRaw).getId().stream().map(e -> ((AIdentifierId)e).getValue().getText()).collect(Collectors.toSet());
  }
  
  private void validate() throws IllegalArgumentException {
    if(!type.contains(LINEAR_OPTION_ID) ^ type.contains(NONLINEAR_OPTION_ID)) {//should be only one
      throw new IllegalArgumentException("nonlinear or linear not specified or both specified in the kernel type");
    }
  }
  
  //use by the solverFactory
  public boolean isKernelType(String id) {
    return optimization.contains(id);
  }
  
  public boolean isLinear() throws IllegalArgumentException {
    return type.contains(LINEAR_OPTION_ID);
  }
  
  public boolean useFlux() {
    return terms.contains(FLUX_OPTION_ID);
  }
  
  public boolean useSource() {
    return terms.contains(SOURCE_OPTION_ID);
  }
  
  public boolean useNCP() {
    return terms.contains(NCP_OPTION_ID);
  }
  
  public boolean usePointSource() {
    return terms.contains(POINTSOURCE_OPTION_ID);
  }
  
  public boolean noTimeAveraging() {
    return optimization.contains(NO_TIME_AVG_OPTION_ID);
  }
  
  //(type: [...], terms: [...], opt: [...])
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("(type: [");
    for(String s : type) {
      sb.append(s);
      sb.append(", ");
    }
    sb.deleteCharAt(sb.length()-2);
    sb.append("], terms: [");
    for(String s : terms) {
      sb.append(s);
      sb.append(", ");
    }
    sb.deleteCharAt(sb.length()-2);
    sb.append("], opt: [");
    for(String s : optimization) {
      sb.append(s);
      sb.append(", ");
    }
    sb.deleteCharAt(sb.length()-2);
    sb.append("])");
    
    return sb.toString();
  }
  
}