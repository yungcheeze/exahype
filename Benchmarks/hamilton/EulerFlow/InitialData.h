#ifndef __INITIAL_DATA_EULERFLOW__
#define __INITIAL_DATA_EULERFLOW__

#include <cstdlib>
#include <iostream>
#include <string>

typedef void (*InitialDataHandler)(const double* const x, double t, double* Q);

extern InitialDataHandler idfunc;

void InitialData(const double* const x, double t, double* Q);

void Polynomial(const double* const x, double t, double* Q);

void ShuVortex2D(const double* const x, double t, double* Q);

void MovingGauss(const double* const x, double t, double* Q);

void DiffusingGauss(const double* const x, double t, double* Q);

void SmoothedSodShockTube(const double* const x, double t, double* Q);

#endif /* __INITIAL_DATA_EULERFLOW__ */
