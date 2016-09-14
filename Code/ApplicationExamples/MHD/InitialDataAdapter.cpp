/**
 * This file is a pure bugfix to circumvent C to FORTRAN parameter passing (especially
 * string passing). As I have not found a better way to push the outcome of
 * constants->getValueAsString() to FORTRAN, the FORTRAN initial data function calls
 * this C function which again in turn calls the FORTRAN functions. This is an ugly
 * workaround for a better future solution.
 **/

#include "MHDSolver.h"
#include "tarch/logging/Log.h"
#include <map> // was an idea, std::map<string, idfunc>
#include <algorithm>

static const char* specfile_initialdata_key = "initialdata";

// I am lazy
#define constants MHDSolver::MHDSolver::constants

typedef void (*idfunc)(double* x, double* Q);

static tarch::logging::Log _log("MHDSolver");

extern "C" {

// FORTRAN functions called by C
void initialalfenwave_(double* x, double* Q);
void initialblast_(double* x, double* Q);
void initialorsagtang_(double* x, double* Q);
void initialrotor_(double* x, double* Q);
	
// C functions called by FORTRAN
void initialdatabyexahypespecfile(double* x, double* Q) {
	if(!constants) {
		logError("InitialDatabyExahyPESpecFile()", "Parser instance is Null");
	} else if(!constants->isValueValidString(specfile_initialdata_key)) {
		logError("InitialDatabyExahyPESpecFile()", "Initial Data Value is not valid");
	} else {
		std::string id_name = constants->getValueAsString(specfile_initialdata_key);
		// c++ string to lower
		std::transform(id_name.begin(), id_name.end(), id_name.begin(), ::tolower);
		
		     if(id_name == "alfenwave") initialalfenwave_(x,Q);
		else if(id_name == "blast")     initialblast_(x,Q);
		else if(id_name == "orsagtang") initialorsagtang_(x,Q);
		else if(id_name == "rotor")     initialrotor_(x,Q);
		else {
			logError("InitialDatabyExahyPESpecFile()", "Unknown Initial Data value");
		}
		return; // success. or so.
	}
	// failure: Do fallback initial data
	logWarning( "InitialDatabyExahyPESpecFile()", "Falling back to AlfenWave");
	initialalfenwave_(x,Q);
}

}/* extern "C" */
