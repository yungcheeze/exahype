// If HDF5 support is not enabled, there will be no implementation of CarpetHDF5Writer
// in the compiled binary.
#ifdef HDF5

#include "exahype/plotters/CarpetHDF5/CarpetHDF5Writer.h"
#include "peano/utils/Loop.h" // dfor
#include <sstream>

using namespace H5;

// my small C++11 to_string-independent workaround.
template <typename T> std::string toString( T Number ) {
	std::ostringstream ss; ss << Number; return ss.str();
}

exahype::plotters::CarpetHDF5Writer::CarpetHDF5Writer(
	const std::string& _filename,
	int _basisSize,
	int _solverUnknowns,
	int _writtenUnknowns,
	const std::string& _select,
	char** _writtenQuantitiesNames,
	bool _oneFilePerTimestep,
	bool _allUnknownsInOneFile)
	:
	_log("ADERDG2CarpetHDF5Writer"),
	solverUnknowns(_solverUnknowns),
	writtenUnknowns(_writtenUnknowns),
	basisFilename(_filename),
	basisSize(_basisSize),
	select(_select),
	component(-100),
	iteration(0),
	writtenQuantitiesNames(_writtenQuantitiesNames),
	#if DIMENSIONS == 2
	// we can also use pointers if you want to add logic here or cannot use
	// preprocessor directives, ie.
	// writtenCellIdx = new kernels::index(basisSize, basisSize, writtenUnknowns) ...etc...
	writtenCellIdx(basisSize, basisSize, writtenUnknowns),
	singleFieldIdx(basisSize, basisSize),
	#else
	writtenCellIdx(basisSize, basisSize, basisSize, writtenUnknowns),
	singleFieldIdx(basisSize, basisSize, basisSize),
	#endif
	single_file(nullptr),
	seperate_files(nullptr),
	oneFilePerTimestep(_oneFilePerTimestep),
	allUnknownsInOneFile(_allUnknownsInOneFile)
	{

	// make sure there are reasonable names everywhere
	for(int u=0; u<writtenUnknowns; u++) {
		char* field_name = writtenQuantitiesNames[u];
		if(!field_name) {
			std::string* replacement_name = new std::string("Q_");
			*replacement_name += toString(u);
			writtenQuantitiesNames[u] = const_cast<char*>(replacement_name->c_str());
		}
	}
	
	
	// this is the dataspace describing how to write a patch/cell/component
	const int dims_rank = DIMENSIONS;
	hsize_t dims[dims_rank];
	std::fill_n(dims, DIMENSIONS, basisSize);
	patch_space = DataSpace(dims_rank, dims);
	
	// this is just a vector of rank 1 with DIMENSIONS entries
	const int tupleDim_rank = 1;
	hsize_t tupleDim_len[] = {DIMENSIONS};
	dtuple = DataSpace(tupleDim_rank, tupleDim_len);
	
	if(!allUnknownsInOneFile) {
		seperate_files = new H5::H5File*[writtenUnknowns];
		std::fill_n(seperate_files, writtenUnknowns, nullptr);
	}

	// open file(s) initially, if neccessary
	if(!oneFilePerTimestep) {
		openH5(0);
	}
}

void exahype::plotters::CarpetHDF5Writer::writeBasicGroup(H5::H5File* file) {
	Group* parameters = new Group(file->createGroup( "/Parameters and Global Attributes" ));
	
	int ranks = 1; // TODO: extend for MPI.
	
	Attribute nioprocs = parameters->createAttribute("nioprocs", PredType::NATIVE_INT, H5S_SCALAR);
	nioprocs.write(PredType::NATIVE_INT, &ranks);
	
	// Todo: add information how much files are printed (oneFilePerTimestep)
}
  
/**
 * We use a H5File** here in order to modify a passed pointer H5File* to a single
 * H5File.
 **/
void exahype::plotters::CarpetHDF5Writer::openH5File(H5::H5File** file, std::string& filename) {
	if(*file) {
		(*file)->close();
		delete *file;
		*file = nullptr;
	}
	
	logInfo("openH5File", "Opening File '"<< filename << "'");
	*file = new H5File( H5std_string(filename.c_str()), H5F_ACC_TRUNC  );
	(*file)->setComment("Created by ExaHyPE");
}
  
/**
 * Opens or switchs the currently active H5 file or the list of H5 files.
 * 
 **/
void exahype::plotters::CarpetHDF5Writer::openH5(int iteration) {
	std::string local_filename, suffix, prefix, sep("-");
	prefix = basisFilename;
	suffix = (oneFilePerTimestep?(sep + "it" + toString(iteration)):"") + ".h5";
	if(allUnknownsInOneFile) {
		local_filename = prefix + suffix;
		openH5File(&single_file, local_filename);
	} else {
		for(int u=0; u < writtenUnknowns; u++) {
			char* writtenQuantityName = writtenQuantitiesNames[u];
			local_filename = prefix + sep + writtenQuantityName + suffix;
			openH5File(&seperate_files[u], local_filename);
		}
	}
}
  
void exahype::plotters::CarpetHDF5Writer::flushH5File(H5::H5File* file) {
	if(file) {
		file->flush(H5F_SCOPE_GLOBAL);
	}
}

void exahype::plotters::CarpetHDF5Writer::startPlotting(double time) {
	component = 0; // CarpetHDF5 wants the components start with 0.
	if(oneFilePerTimestep) {
		openH5(iteration);
	}
}
  
void exahype::plotters::CarpetHDF5Writer::finishPlotting() {
	// flush in any case, even when closing the file. Cf. http://stackoverflow.com/a/31301117
	if(allUnknownsInOneFile) {
		flushH5File(single_file);
	} else {
		for(int u=0; u < writtenUnknowns; u++) {
			flushH5File(seperate_files[u]);
		}
	}
	iteration++;
}

/**
 * This is 2D and 3D, allows several unknowns, named fields and all that.
 **/
void exahype::plotters::CarpetHDF5Writer::plotPatch(
      const dvec& offsetOfPatch, const dvec& sizeOfPatch, const dvec& dx,
      double* mappedCell, double timeStamp) {
	for(int writtenUnknown=0; writtenUnknown < writtenUnknowns; writtenUnknown++) {
		H5::H5File* target = allUnknownsInOneFile ? single_file : seperate_files[writtenUnknown];
		plotPatchForSingleUnknown(offsetOfPatch, sizeOfPatch, dx, mappedCell, timeStamp, writtenUnknown, target);
	} // for writtenUnknown
	component++;
}

void exahype::plotters::CarpetHDF5Writer::plotPatchForSingleUnknown(
      const dvec& offsetOfPatch, const dvec& sizeOfPatch, const dvec& dx,
      double* mappedCell, double timeStamp,
      int writtenUnknown, H5::H5File* target) {
	// in CarpetHDF5, the field name *must* contain a "::"
	std::string name("ExaHyPE::");
	char* field_name = writtenQuantitiesNames[writtenUnknown];
	name += field_name ? field_name : "miserable-failure";
	
	char component_name[100];
	sprintf(component_name, "%s it=%d tl=0 m=0 rl=0 c=%d", name.c_str(), iteration, component);
	
	// 1) Compose a continous storage which is suitable.
	// 2) (Probably) Transpose the data.
	// TODO: I'm sure HDF5 provides a more performant way to interpret the different data layout.
	double *componentPatch = new double[singleFieldIdx.size];
	dfor(i,basisSize) {
		#if DIMENSIONS==2
		componentPatch[singleFieldIdx(i(1),i(0))] = mappedCell[writtenCellIdx(i(0),i(1),writtenUnknown)];
		#else
		componentPatch[singleFieldIdx(i(2),i(1),i(0))] = mappedCell[writtenCellIdx(i(0),i(1),i(2),writtenUnknown)];
		#endif
	}
	
	DataSet table = target->createDataSet(component_name, PredType::NATIVE_FLOAT, patch_space);
	table.write(componentPatch, PredType::NATIVE_DOUBLE);
	
	delete[] componentPatch;
	
	// write all meta information about the table

	Attribute origin = table.createAttribute("origin", PredType::NATIVE_DOUBLE, dtuple);
	origin.write(PredType::NATIVE_DOUBLE, offsetOfPatch.data());
	
	dvec iorigin_data = 0.0;
	Attribute iorigin = table.createAttribute("iorigin", PredType::NATIVE_INT, dtuple);
	iorigin.write(PredType::NATIVE_INT, iorigin_data.data());
	
	int level_data = 0; // TODO: Read out real cell level
	Attribute level = table.createAttribute("level", PredType::NATIVE_INT, H5S_SCALAR);
	level.write(PredType::NATIVE_INT, &level_data);
	
	Attribute timestep = table.createAttribute("timestep", PredType::NATIVE_INT, H5S_SCALAR);
	timestep.write(PredType::NATIVE_INT, &iteration);

	Attribute time = table.createAttribute("time", PredType::NATIVE_DOUBLE, H5S_SCALAR);
	time.write(PredType::NATIVE_DOUBLE, &timeStamp);

	// dx in terms of Cactus: Real seperation from each value
	Attribute delta = table.createAttribute("delta", PredType::NATIVE_FLOAT, dtuple);
	delta.write(PredType::NATIVE_DOUBLE, dx.data()); // issue: conversion from double to float
	
	const int max_string_length = 60;
	StrType t_str = H5::StrType(H5::PredType::C_S1, max_string_length); // todo: use name.size().
	Attribute aname = table.createAttribute("name", t_str, H5S_SCALAR);
	aname.write(t_str, name.c_str());
}


#endif /* HDF5 */
