package eu.exahype.solvers;

import java.io.BufferedWriter;
import java.io.IOException;
import java.lang.IllegalArgumentException;
import java.util.Set;
import java.util.List;

// template engine
import minitemp.Context;
import minitemp.TemplateEngine;

import eu.exahype.CodeGeneratorHelper;
import eu.exahype.io.IOUtils;
import eu.exahype.io.SourceTemplate;



public class OptimisedADERDG implements Solver {
  
  public static final String Identifier = "optimised"; //expect the spec file to have the form "optimised::options"

  //Internal states
  //--------------- 
  private String _solverName;
  private Context context;
  private TemplateEngine templateEngine;
  
  private boolean _isValid; //if false the solverFactory will return null

  public OptimisedADERDG(String projectName, String solverName, int dimensions, int numberOfVariables, int numberOfParameters, Set<String> namingSchemeNames,
      int order,String microarchitecture, boolean enableProfiler, boolean enableDeepProfiler, boolean hasConstants, List<String> options) 
      throws IOException, IllegalArgumentException {    
    
    _solverName         = solverName;
    
    if(!options.contains(LINEAR_OPTION_ID) ^ options.contains(NONLINEAR_OPTION_ID)) {//should be only one
      throw new IllegalArgumentException("nonlinear or linear not specified or both specified in the options ("+options+")");
    }
    final boolean isLinear           = options.contains(LINEAR_OPTION_ID);
    final boolean useFlux            = options.contains(FLUX_OPTION_ID);
    final boolean useSource          = options.contains(SOURCE_OPTION_ID);
    final boolean useNCP             = options.contains(NCP_OPTION_ID);
    final boolean noTimeAveraging    = options.contains(NO_TIME_AVG_OPTION_ID); 
    
    templateEngine = new TemplateEngine();
    context = new Context();
    
    //String
    context.put("project"           , projectName);
    context.put("solver"            , solverName);
    context.put("abstractSolver"    , getAbstractSolverName());
    
    //int
    context.put("dimensions"        , dimensions);
    context.put("order"             , order);
    context.put("numberOfVariables" , numberOfVariables);
    context.put("numberOfParameters", numberOfParameters);
    
    //boolean
    context.put("enableProfiler"    , enableProfiler);
    context.put("enableDeepProfiler", enableDeepProfiler);
    context.put("hasConstants"      , hasConstants);
    context.put("isLinear"          , isLinear);
    context.put("useFlux"           , useFlux);
    context.put("useSource"         , useSource);
    context.put("useNCP"            , useNCP);
    context.put("noTimeAveraging"   , noTimeAveraging);
    
    //boolean as String
    context.put("useFlux_s"         , boolToTemplate(useFlux));
    context.put("useSource_s"       , boolToTemplate(useSource));
    context.put("useNCP_s"          , boolToTemplate(useNCP));

    //generate the optimised kernel, can throw IOException
    final String optKernelPath = CodeGeneratorHelper.getInstance().invokeCodeGenerator(projectName, solverName, numberOfVariables, numberOfParameters, order, isLinear, dimensions,
        microarchitecture, enableDeepProfiler, useFlux, useSource, useNCP, noTimeAveraging);
    final String optNamespace = CodeGeneratorHelper.getInstance().getNamespace(projectName, solverName);
    
    context.put("optKernelPath"     , optKernelPath);
    context.put("optNamespace"      , optNamespace);
    
    
    //TODO JM support for loop (for writeAbstractImplementation)
    String namingSchemes = "";
    for (String name : namingSchemeNames) {
      namingSchemes += "    " + "class "+name.substring(0, 1).toUpperCase() + name.substring(1) + ";\n";
    }
    context.put("namingSchemes"     , namingSchemes);
    
    
    //TODO JMG: linear kernels unsupported for now
    if(isLinear) 
      throw new IllegalArgumentException("Linear kernels not supported yet");
  }
    
  @Override
  public String getSolverName() {
    return _solverName;
  }
  
  private String getAbstractSolverName() {
    return "Abstract"+getSolverName();
  }
  
  private String boolToTemplate(boolean b) {
    return b? "true" : "false";
  }
  
  @Override
  public void writeHeader(java.io.BufferedWriter writer)
      throws IOException, IllegalArgumentException {
	  final String template = IOUtils.convertRessourceContentToString("eu/exahype/solvers/templates/OptimisedADERDGSolverHeader.template");
	  writer.write(templateEngine.render(template, context));
  }

  
  @Override
  public void writeAbstractHeader(java.io.BufferedWriter writer)
      throws java.io.IOException, IllegalArgumentException {      
    final String template = IOUtils.convertRessourceContentToString("eu/exahype/solvers/templates/AbstractOptimisedADERDGSolverHeader.template");
    writer.write(templateEngine.render(template, context));
  }
  
  @Override
  public void writeAbstractImplementation(java.io.BufferedWriter writer) throws java.io.IOException, IllegalArgumentException {

    final String template = IOUtils.convertRessourceContentToString("eu/exahype/solvers/templates/AbstractOptimisedADERDGSolverImplementation.template"); 
    writer.write(templateEngine.render(template, context));
  }
  
  //TODO JMG move to template engine
  @Override
  public void writeUserImplementation(java.io.BufferedWriter writer) throws java.io.IOException, IllegalArgumentException {
/*     SourceTemplate content = SourceTemplate.fromRessourceContent(
        "eu/exahype/solvers/templates/GenericADERDGSolverInCUserCode.template");
    
    content.put("Project", _projectName);
    content.put("Solver", _solverName);
    
    content.put("Elements",  String.valueOf( _numberOfParameters+_numberOfVariables));
    content.put("Dimensions",String.valueOf(_dimensions));

    String SolverInitSignatureExtension = "";
    if (_hasConstants) {
        SolverInitSignatureExtension = ", exahype::Parser::ParserView& constants";
    }
    content.put("SolverInitSignatureExtension", SolverInitSignatureExtension);
    //
    String solverConstructorArgumentExtension  = "";
    String solverConstructorSignatureExtension = "";
    String SolverInitCallExtension             = "";
    if (_enableProfiler) {
      solverConstructorSignatureExtension += ", std::unique_ptr<exahype::profilers::Profiler> profiler";
      solverConstructorArgumentExtension  += ", std::move(profiler)";
    }
    if (_hasConstants) {
      solverConstructorSignatureExtension += ", exahype::Parser::ParserView& constants";
       SolverInitCallExtension = ", constants";
    }

    content.put("SolverInitCallExtension",SolverInitCallExtension);
    content.put("SolverConstructorSignatureExtension", solverConstructorSignatureExtension);
    content.put("SolverConstructorArgumentExtension", solverConstructorArgumentExtension);
    
    // user functions
    int digits = String.valueOf(_numberOfVariables + _numberOfParameters).length();
    
    String adjustSolution = "  // State variables:\n";
    for (int i = 0; i < _numberOfVariables; i++) {
      adjustSolution += "  Q[" + String.format("%" + digits + "d", i) + "] = 0.0;";
      if (i<_numberOfVariables-1) adjustSolution += "\n";
    }
    if (_numberOfParameters>0) {
      adjustSolution += "  // Material parameters:\n";
      for (int i = 0; i < _numberOfParameters; i++) {
        adjustSolution += "  Q[" + String.format("%" + digits + "d", _numberOfVariables+i) + "] = 0.0;";
        if (i<_numberOfParameters-1) adjustSolution += "\n";
      }
    }

    String eigenvalues = "";
    for (int i = 0; i < _numberOfVariables; i++) {
      eigenvalues += "  lambda[" + String.format("%" + digits + "d", i) + "] = 0.0;";
      if (i<_numberOfVariables-1) eigenvalues += "\n";
    }

    String flux = "";
    for (int d=0; d<_dimensions; ++d) {
      for (int i = 0; i < _numberOfVariables; i++) {
        flux += "  F["+d+"][" + String.format("%" + digits + "d", i) + "] = 0.0;";
        if (i<_numberOfVariables-1) flux += "\n";
      }
      if (d<_dimensions-1) {
        flux += "\n\n";    
      }
    }
    
    String source = "";
    for (int i = 0; i < _numberOfVariables; i++) {
      source += "  S[" + String.format("%" + digits + "d", i) + "] = 0.0;";
      if (i<_numberOfVariables-1) source += "\n";
    }
    
    String boundaryValues = "";
    for (int i = 0; i < _numberOfVariables; i++) {
      boundaryValues += "  stateOut[" + String.format("%" + digits + "d", i) + "] = 0.0;";
      if (i<_numberOfVariables-1) boundaryValues += "\n";
    }
    boundaryValues += "\n\n";
    for (int i = 0; i < _numberOfVariables; i++) {
      boundaryValues += "  fluxOut[" + String.format("%" + digits + "d", i) + "] = 0.0;";
      if (i<_numberOfVariables-1) boundaryValues += "\n";
    }
    
    String ncp = "";
    for (int i = 0; i < _numberOfVariables; i++) {
      ncp += "  BgradQ[" + String.format("%" + digits + "d", i) + "] = 0.0;";
      if (i<_numberOfVariables-1) ncp += "\n";
    }
    
    String matrixb = "";
    for (int i = 0; i < _numberOfVariables*_numberOfVariables; i++) {
      matrixb += "  Bn[" + String.format("%" + digits + "d", i) + "] = 0.0;";
      if (i<_numberOfVariables*_numberOfVariables-1) matrixb += "\n";
    }
    
    content.put("AdjustedSolutionValues",adjustSolution);
    content.put("Eigenvalues",eigenvalues);
    content.put("Flux",flux);
    content.put("Source",source);
    content.put("BoundaryValues",boundaryValues);
    content.put("NonConservativeProduct",ncp);
    content.put("MatrixB",matrixb);
    
    writer.write(content.toString()); */
  }
  
  @Override
  public boolean supportsVariables() {
    return true;
  }
}
