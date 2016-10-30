# TwoPunctures

TwoPunctures is an Initial Data code for General Relativity user applications. It creates
consistent initial data for Black Hole spacetimes.

This code was extracted from EinsteinToolkit and is currently hosted at

  https://bitbucket.org/relastro/twopunctures-standalone

As it is a very tiny code, it can probably just imported to ExaHyPE some day.

# Why this place

This code cannot live inside the Applications directory as the full repository contains
examples including a `main()`. Instead, we build the static library manually and link
it to the ExaHyPE binary.

# Dependencies

TwoPunctures requires the GNU Scientific Library (GSL) to compile and execute.
The GSL itself needs BLAS. You end up with a linking like `-lgsl -lgslcblas -lm`

