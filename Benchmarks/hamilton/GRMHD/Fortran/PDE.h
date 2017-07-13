#ifndef __EXAHYPE_USER_PDE__
#define __EXAHYPE_USER_PDE__

// Fortran functions:
extern "C" {
void minimumtreedepth_(int* depth);
void hastoadjustsolution_(double* t, bool* refine);
void adjustedsolutionvalues_(const double* const x,const double* w,const double* t,const double* dt,double* Q);
void pdeflux_(double* Fx, double* Fy, double* Fz, const double* const Q);
void pdesource_(double* S, const double* const Q);
void pdeeigenvalues_(double* lambda, const double* const Q, double* nv);
void registerinitialdata_(const char* const id_name, int* id_name_len);

void pdencp_(double* BgradQ, const double* const Q, const double* const gradQ);
void pdematrixb_(double* Bn, const double* const Q, double* nv);

// C2P-GRMHD.f90
void pdeprim2cons_(double* /* OUT */ Q, double* /* IN */ V);
void pdecons2prim_(double* /* OUT */ V, double* /* IN */ Q, int* iErr);

// InitialDataAlfenWave.f90
void alfenwave_(const double* x,const double* const t,  double* Q);


}/* extern "C" */

#endif /* __EXAHYPE_USER_PDE__ */
