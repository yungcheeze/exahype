#ifndef __INITIAL_DATA_EULERFLOW__
#define __INITIAL_DATA_EULERFLOW__

#include <iostream>
#include <cstdlib>
#include <string>

/* This is where we get our parameter decision from: */
inline const std::string& getEnvParameter(const char* key) {
	if(char* env = std::getenv(key))
		return env;
	else return "";
}

void InitialData(const double* const  x, double* Q, double t);

void ShuVortex2D(const double* const  x, double* Q, double t);

#endif /* __INITIAL_DATA_EULERFLOW__ */