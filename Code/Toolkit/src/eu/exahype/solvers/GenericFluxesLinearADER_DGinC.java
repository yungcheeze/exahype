package eu.exahype.solvers;

public final class GenericFluxesLinearADER_DGinC extends GenericFluxesADER_DGinC {
  public static final String Identifier = "generic::fluxes::linear";

  public GenericFluxesLinearADER_DGinC(int dimensions, int numberOfUnknowns, int numberOfParameters,
      int order, boolean enableProfiler, boolean hasConstants) {
    super(dimensions, numberOfUnknowns, numberOfParameters, order, enableProfiler, hasConstants);
  }

  @Override
  public final boolean isLinear() {
    return true;
  }

  @Override
  public final boolean isFortran() {
    return false;
  }
}
