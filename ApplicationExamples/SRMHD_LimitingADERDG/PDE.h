#ifndef __EXAHYPE_USER_PDE__
#define __EXAHYPE_USER_PDE__

// Fortran functions:
extern "C" {
void minimumtreedepth_(int* depth);
void hastoadjustsolution_(double* t, bool* refine);
void adjustedsolutionvalues_(const double* const x,const double* w,const double* t,const double* dt,double* Q);
void pdeflux_(double* F, const double* const Q);
void pdesource_(double* S, const double* const Q);
void pdeeigenvalues_(double* lambda, const double* const Q, double* nv);
void registerinitialdata_(const char* const id_name, int* id_name_len);

void pdeboundaryvalues_(
	const double* const x, const double* const t, const double* const dt,
	const int* const faceIndex, double* nv,
	const double* const fluxIn, const double* const stateIn,
	double* fluxOut, double* stateOut);

}/* extern "C" */

#endif /* __EXAHYPE_USER_PDE__ */
