#ifndef __MASSACCRETIONRATE_CPP__
#define __MASSACCRETIONRATE_CPP__

/**
 * C definitions to use the GRMHD functions
 **/

extern "C" {

// FORTRAN functions called by C
void massaccretionrate_(double* /* IN */ Q, double* /* OUT */ masschange);

} /* extern "C" */
#endif /* __MASSACCRETIONRATE_CPP__ */
