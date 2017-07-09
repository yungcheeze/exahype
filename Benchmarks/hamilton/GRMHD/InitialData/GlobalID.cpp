#include "InitialData/InitialData.h"

idobj *id = nullptr; // storage

// a function in InitialData.cpp which prepares the ID, both accessible
// from a pure ADERDG, pure FV or limiting application
bool prepare_id(std::string idname) {
	if(id) {
		// id already prepared
		return true;
	}
	
	if(idname == "AlfenWave") {
		id = new fortranid();
		return true;
	} else if(idname == "PizzaTOV") {
		id = new pizzatov();
		return true;
	} else if(idname == "RNSID") {
		id = new rnsid();
		return true;
	}
	return false; // no success
}
