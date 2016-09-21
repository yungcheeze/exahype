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

#ifndef _EXAHYPE_TEST_OPT_KERNELS_H_
#define _EXAHYPE_TEST_OPT_KERNELS_H_

#include "matrices.h"

// for order 3..8

void picardLoop3(double* restrict lqh, double* restrict lFh, const double* restrict const luh, const double* dx, const double dt);
void picardLoop4(double* restrict lqh, double* restrict lFh, const double* restrict const luh, const double* dx, const double dt);
void picardLoop5(double* restrict lqh, double* restrict lFh, const double* restrict const luh, const double* dx, const double dt);
void picardLoop6(double* restrict lqh, double* restrict lFh, const double* restrict const luh, const double* dx, const double dt);
void picardLoop7(double* restrict lqh, double* restrict lFh, const double* restrict const luh, const double* dx, const double dt);
void picardLoop8(double* restrict lqh, double* restrict lFh, const double* restrict const luh, const double* dx, const double dt);


#endif