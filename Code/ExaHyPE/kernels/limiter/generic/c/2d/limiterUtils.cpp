/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/

#include "../../Limiter.h"
 
#if DIMENSIONS == 2
 
namespace kernels {
namespace limiter {
namespace generic {
namespace c {

//Fortran (Limiter.f90): GetLobattoData
//TODO JMG
double* getGaussLobattoData(const double* const luh, const int numberOfVariables, const int basisSize) {
  double* lob = new double[basisSize*basisSize*numberOfVariables]; //Fortran ref: lob(nVar,nDOF(1),nDOF(2),nDOF(3))
  
  return lob;
}

//Fortran (Limiter.f90): GetSubcellData
//TODO JMG
double* getFVMData(const double* const luh, const int numberOfVariables, const int basisSize, int& basisSizeLim) {
  
  basisSizeLim = getLimBasisSize(basisSize);
 
  double* lim = new double[basisSizeLim*basisSizeLim*numberOfVariables]; //Fortran ref: lim(nVar,nSubLimV(1),nSubLimV(2),nSubLimV(3))

  return lim;
}


} // namespace c
} // namespace generic
} // namespace limiter
} // namespace kernel

#endif //DIMENSIONS == 2