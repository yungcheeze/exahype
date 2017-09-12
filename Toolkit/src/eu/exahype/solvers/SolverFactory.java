package eu.exahype.solvers;

import java.util.Set;
import java.util.List;
import java.util.Arrays;

public class SolverFactory {
  private String _projectName;
  private int _dimensions;
  private boolean _enableProfiler;
  private boolean _enableDeepProfiler;
  private String _microarchitecture;
  private String _pathToLibxsmm;

  
  public SolverFactory(
      String projectName,
      int dimensions,
      boolean enableProfiler,
      boolean enableDeepProfiler,
      String microarchitecture) {
    _projectName = projectName;  
    _dimensions = dimensions;
    _enableProfiler = enableProfiler;
    _enableDeepProfiler = enableDeepProfiler;
    _microarchitecture = microarchitecture;
  }
  
  public Solver createADERDGSolver(String solvername, String kernel,boolean isFortran,int numberOfVariables,int numberOfParameters,Set<String> namingSchemeNames,int order,boolean hasConstants) {
    String generalKernel = kernel.substring(0, kernel.lastIndexOf("::"));
    boolean isLinear     = kernel.substring(kernel.lastIndexOf("::")).equalsIgnoreCase("::linear");
    
    if (isFortran && kernel.equals( eu.exahype.solvers.UserDefinedADER_DGinFortran.Identifier )) {
      return new eu.exahype.solvers.UserDefinedADER_DGinFortran(_projectName, solvername);
    }
    else if (!isFortran && kernel.equals( eu.exahype.solvers.UserDefinedADER_DGinC.Identifier )) {
      return new eu.exahype.solvers.UserDefinedADER_DGinC(_projectName, solvername, numberOfVariables,
          numberOfParameters, order, hasConstants, _enableProfiler);
    }
    else if (generalKernel.equals( eu.exahype.solvers.GenericADERDG.Identifier )) {
      return new eu.exahype.solvers.GenericADERDG(_projectName, solvername, _dimensions,
          numberOfVariables, numberOfParameters, namingSchemeNames, order, _enableProfiler, hasConstants, isLinear, isFortran );
    }
    // TODO JMG Clean
    // else if (!isFortran && kernel.equals( eu.exahype.solvers.OptimisedFluxesLinearADER_DGinC.Identifier )) {
      // return new eu.exahype.solvers.OptimisedFluxesLinearADER_DGinC(_dimensions,
          // numberOfVariables, numberOfParameters, order, _microarchitecture,
          // _enableProfiler, hasConstants);
    // }
    else if (!isFortran && generalKernel.startsWith( eu.exahype.solvers.OptimisedADERDG.Identifier )) {
      eu.exahype.solvers.OptimisedADERDG solver =  new eu.exahype.solvers.OptimisedADERDG(_projectName, solvername, _dimensions,
          numberOfVariables, numberOfParameters, namingSchemeNames, order, _microarchitecture,
          _enableProfiler, _enableDeepProfiler, hasConstants, Arrays.asList(kernel.split("::")));
      return (solver.isValid() ? solver : null);
    }
    else if (!isFortran && kernel.equals( eu.exahype.solvers.KernelEuler2d.Identifier )) {
      return new eu.exahype.solvers.KernelEuler2d(_projectName, solvername);
    }
    
    return null;
  }
  
  public Solver createFiniteVolumesSolver(String solvername, String kernel,boolean isFortran,int numberOfVariables,int numberOfParameters,Set<String> namingSchemeNames,int patchSize,boolean hasConstants) {
    if (isFortran && kernel.equals( eu.exahype.solvers.UserDefinedFiniteVolumesinFortran.Identifier )) {
      return new eu.exahype.solvers.UserDefinedFiniteVolumesinFortran(_projectName, solvername,_dimensions,numberOfVariables, numberOfParameters, patchSize, _enableProfiler, hasConstants);
    }
    if (!isFortran && kernel.equals( eu.exahype.solvers.UserDefinedFiniteVolumesinC.Identifier )) {
      return new eu.exahype.solvers.UserDefinedFiniteVolumesinC(_projectName, solvername,_dimensions,numberOfVariables, numberOfParameters, patchSize, _enableProfiler, hasConstants);
    }
    if (!isFortran && kernel.equals( eu.exahype.solvers.GenericFiniteVolumesGodunovInC.Identifier )) {
      return new eu.exahype.solvers.GenericFiniteVolumesGodunovInC(_projectName, solvername,_dimensions,numberOfVariables, numberOfParameters, namingSchemeNames, patchSize, _enableProfiler, hasConstants);
    }
    if (!isFortran && kernel.equals( eu.exahype.solvers.GenericFiniteVolumesMUSCLHancockInC.Identifier )) {
        return new eu.exahype.solvers.GenericFiniteVolumesMUSCLHancockInC(_projectName, solvername,_dimensions,numberOfVariables, numberOfParameters, namingSchemeNames, patchSize, _enableProfiler, hasConstants);
    }

    return null;
  }
}
