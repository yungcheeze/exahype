/**
 * This are the function signatures for the actual FORTRAN kernels,
 * masqueraded by the C++ functions.
 * 
 * 
 **/

#define nVar 3
#define d    3

typedef double*  fdvec; ///< Fortran double vector
typedef int      fsize; ///< Fortran vector length information
typedef int*     fint;  ///< Fortran integer

PDEFlux_(fdvec F, fdvec Q, fsize fLength, fsize qLength);
PDEEigenvalues_(fdvec Lambda, fdvec Q, fdvec nv, fsize LambdaLength, fsize QLength, fsize nvLength);
PDECons2Prim_(fdvec V, fdvec Q, fint iErr, fsize Vlen, fsize Qlen);
PDEPrim2Cons_(fdvec Q, fdvec V, fsize Qlen, fsize Vlen)