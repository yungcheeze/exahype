#include "InitialData/InitialData.h"

idobj *id = nullptr; // storage

// a function in InitialData.cpp which prepares the ID, both accessible
// from a pure ADERDG, pure FV or limiting application
void prepare_id() {
	if(id) {
		// id already prepared
		return;
	}
	//id = new fortranid();
	id = new pizzatov();
}
