# Thoughts
[1] The initial values and boundary conditions should move to the
application folder since these are user dependent.

[2] For an optimized  ExaHyPE production run, there will be no 
PDEDescription.cpp/h file specifying fluxes and eigenvalues.
Most/all kernel functions will only contain optimized assembler code.  

[3] To enable rapid prototyping, we should however allow the user to
use some  default kernels and to provide his flux and eigenvalue 
definitions separately.
 
Example:
   use-optimized-kernels    = ${EXAHYPE_DIR}/euler/compressible_euler 
vs.
   use-user-pde-description = ${PROJECT}/PDEDescription.cpp 
   
[4] Mappings, quadrature points etc. should not be part of the 
ExaHyPE engine.


