// If HDF5 support is not enabled, there will be no implementation of CarpetHDF5Writer
// in the compiled binary.
#ifdef HDF5

#include "exahype/plotters/CarpetHDF5/CarpetHDF5Writer.h"
#include "peano/utils/Loop.h" // dfor
#include <sstream>
#include <stdexcept> // why not

typedef tarch::la::Vector<DIMENSIONS, double> dvec;
typedef tarch::la::Vector<DIMENSIONS, bool> boolvec;
typedef tarch::la::Vector<DIMENSIONS, int> ivec;

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
	
	// default values for init(...)-level runtime parameters:
	dim(DIMENSIONS),
	slicer(nullptr),
	component(-100),
	iteration(0),
	writtenQuantitiesNames(_writtenQuantitiesNames),
	single_file(nullptr),
	seperate_files(nullptr),
	oneFilePerTimestep(_oneFilePerTimestep),
	allUnknownsInOneFile(_allUnknownsInOneFile)
	{

	// while CarpetHDF5 files also can be 1D, this writer currently only understands 2D and 3D.
	assert(dim == 2 || dim == 3);
	
	// parse the selection string
	dvec r;
	r(0) = Parser::getValueFromPropertyString(select, "x");
	r(1) = Parser::getValueFromPropertyString(select, "y");
	#if DIMENSIONS==3
	r(2) = Parser::getValueFromPropertyString(select, "z");
	#endif
	ivec ron;
	// NaN means the property was not present in the select string
	for(int i=0; i<DIMENSIONS; i++) { ron(i) = (r(i)!=r(i)) ? 0 : 1; }
	int numberOfRequestedDimensionalReductions = tarch::la::sum(ron);
	int targetDimension = dim - numberOfRequestedDimensionalReductions;
	
	if(targetDimension == 0) {
		logError("SlicingParser", "Running in " << dim << " dimensions, you requested 0-dimensional output slicing '" << select << "' which not supported by the CarpetHDF5Writer. Use Seismograms instead.");
		throw std::invalid_argument("0-dimensional output not supported");
	}
	
	if(targetDimension == 1) {
		logError("SlicingParser", "Running in " << dim << " dimensions, you requested 1-dimensional output slicing '" << select << "', however the CarpetHDF5 plotter does not yet support 1D output.");
		throw std::invalid_argument("1-dimensional output not YET supported");
	}
	
	if(targetDimension != dim) {
		slicer = new exahype::plotters::CartesianSlicer(r, ron);
		dim = targetDimension;
	}

	// now where we now how many dimensions we have to deal with
	writtenCellIdx = new kernels::index(basisSize, basisSize, dim==3 ? basisSize : writtenUnknowns, dim == 3 ? writtenUnknowns : 1);
	singleFieldIdx = new kernels::index(basisSize, basisSize, dim==3 ? basisSize : 1);
	allFieldsSize = writtenCellIdx->size;
	singleFieldSize = singleFieldIdx->size;
	
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
	hsize_t *dims = new hsize_t[dim];
	std::fill_n(dims, dim, basisSize);
	patch_space = DataSpace(dim, dims);
	
	// this is just a vector of rank 1 with dim entries
	const int tupleDim_rank = 1;
	hsize_t tupleDim_len[tupleDim_rank];
	std::fill_n(tupleDim_len, tupleDim_rank, dim);
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
	
	// Important todo: List all fields which go into this file.
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
	writeBasicGroup(*file);
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
	// TODO: I'm sure HDF5 provides a more performant way to interpret the different data layout.
	double *componentPatch = new double[singleFieldIdx->size];
	dfor(i,basisSize) {
		if(dim==2) {
			componentPatch[singleFieldIdx->get(i(1),i(0))] = mappedCell[writtenCellIdx->get(i(1),i(0),writtenUnknown)];
		} else {
			componentPatch[singleFieldIdx->get(i(2),i(1),i(0))] = mappedCell[writtenCellIdx->get(i(2),i(1),i(0),writtenUnknown)];
		}
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
