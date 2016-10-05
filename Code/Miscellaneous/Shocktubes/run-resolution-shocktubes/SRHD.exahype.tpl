/**

 Special Relativistic Hydrodynamics
 Shocktube with different resolutions

 TEMPLATE.

 */
exahype-project  SRHD

  peano-path                 = ./Peano/peano
  tarch-path                 = ./Peano/tarch
  multiscalelinkedcell-path  = ./Peano/multiscalelinkedcell
  sharedmemoryoracles-path   = ./Peano/sharedmemoryoracles
  exahype-path               = ./ExaHyPE
  output-directory           = ./Applications/SRHD
  architecture               = noarch

  computational-domain
    dimension                = 2
    width                    = 1.0, 1.0
    offset                   = 0.0, 0.0
    end-time                 = 0.5
  end computational-domain
  
  shared-memory
    identifier               = dummy
    cores                    = 30
    properties-file          = sharedmemory.properties
  end shared-memory
  
  optimisation
    fuse-algorithmic-steps        = on
    fuse-algorithmic-steps-factor = 0.99
  end optimisation

  solver ADER-DG SRHDSolver
    variables          = 5
    parameters         = 0
    order              = 1
    // @num_points@ grid points = 1/@num_points@ = @dx@
    maximum-mesh-size  = @dxeps@
    time-stepping      = global
    kernel             = generic::fluxes::nonlinear
    language           = C

    plot vtk::binary
      time     = 0.0
      repeat   = 0.004
      output   = ./output-l@level@/solution
      select   = {all}
    end plot
  end solver

end exahype-project
