#include "exahype/plotters/CarpetHDF5/ADERDG2CarpetHDF5.h"

#include <cstdlib>
#include <stdio.h>
#include <sstream>
#include <memory>


std::string exahype::plotters::ADERDG2CarpetHDF5::getIdentifier() {
	return std::string("hdf5::Cartesian::Vertices");
}


// my small C++11 to_string-independent workaround.
template <typename T> std::string toString( T Number ) {
	std::ostringstream ss; ss << Number; return ss.str();
}

typedef tarch::la::Vector<DIMENSIONS, double> dvec;
typedef tarch::la::Vector<DIMENSIONS, int> ivec;

tarch::logging::Log exahype::plotters::ADERDG2CarpetHDF5::_log("exahype::plotters::ADERDG2CarpetHDF5");


#ifndef HDF5
/*************************************************************************************************
 * ADERDG2CarpetHDF5 Dummy implementation in case HDF5 support is skipped.
 *************************************************************************************************/

exahype::plotters::ADERDG2CarpetHDF5::ADERDG2CarpetHDF5(
  exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing) : Device(postProcessing) {
	printf("ERROR: Compile with HDF5, otherwise you cannot use the HDF5 plotter.\n");
	abort();
}

// all other methods are stubs
exahype::plotters::ADERDG2CarpetHDF5::~ADERDG2CarpetHDF5() {}
void exahype::plotters::ADERDG2CarpetHDF5::init(const std::string& filename, int orderPlusOne, int solverUnknowns, int writtenUnknowns, const std::string& select) {}
void exahype::plotters::ADERDG2CarpetHDF5::plotPatch(const int cellDescriptionsIndex, const int element) {}
void exahype::plotters::ADERDG2CarpetHDF5::plotPatch(const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,double timeStamp) {}
void exahype::plotters::ADERDG2CarpetHDF5::startPlotting(double time) {}
void exahype::plotters::ADERDG2CarpetHDF5::finishPlotting() {}

#else
#include "exahype/plotters/CarpetHDF5/CarpetHDF5Writer.h"
#include "kernels/KernelUtils.h" // indexing
#include "peano/utils/Loop.h" // dfor
#include "kernels/DGMatrices.h"
#include "exahype/solvers/ADERDGSolver.h"
#include "kernels/DGBasisFunctions.h"
#include "tarch/logging/Log.h"
#include <sstream>

/*************************************************************************************************
 * ADERDG2CarpetHDF5 non-dummy implementation
 *************************************************************************************************/

exahype::plotters::ADERDG2CarpetHDF5::ADERDG2CarpetHDF5(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing) :
    Device(postProcessing) { writer = nullptr; }

exahype::plotters::ADERDG2CarpetHDF5::~ADERDG2CarpetHDF5() {
	if(writer) delete writer;
}

void exahype::plotters::ADERDG2CarpetHDF5::init(const std::string& filename, int basisSize, int solverUnknowns, int writtenUnknowns, const std::string& select) {
	bool oneFilePerTimestep = true;
	bool allUnknownsInOneFile = true;

	// Determine names of output fields
	char **writtenQuantitiesNames = new char*[writtenUnknowns];
	std::fill_n(writtenQuantitiesNames, writtenUnknowns, nullptr);
	_postProcessing->writtenQuantitiesNames(writtenQuantitiesNames);
	
	writer = new exahype::plotters::CarpetHDF5Writer(filename, basisSize, solverUnknowns, writtenUnknowns, select,
		writtenQuantitiesNames, oneFilePerTimestep, allUnknownsInOneFile);	

	if(writer->slicer) {
		logInfo("init", "Plotting selection "<<writer->slicer->toString()<<" to Files "<<filename);
	}
}

void exahype::plotters::ADERDG2CarpetHDF5::plotPatch(
        const int cellDescriptionsIndex,
        const int element) {
  auto& aderdgCellDescription = exahype::solvers::ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

  if (aderdgCellDescription.getType()==exahype::solvers::ADERDGSolver::CellDescription::Type::Cell) {
    double* solverSolution = DataHeap::getInstance().getData(aderdgCellDescription.getSolution()).data();

    plotPatch(
        aderdgCellDescription.getOffset(),
        aderdgCellDescription.getSize(), solverSolution,
        aderdgCellDescription.getCorrectorTimeStamp());
  }
}

void exahype::plotters::ADERDG2CarpetHDF5::plotPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    double* u,
    double timeStamp) {

    if(writer->slicer && !writer->slicer->isPatchActive(offsetOfPatch, sizeOfPatch)) {
	return;
    }

    const int basisSize = writer->basisSize;
    const int order     = basisSize - 1;
    const dvec dx       = 1./order * sizeOfPatch;

    // TODO: if we knew that plotting would be serial, we could move *mappedCell to a class property.
    double* mappedCell  = new double[writer->allFieldsSize];

    if(writer->slicer) {    
	interpolateSlicedPatch(offsetOfPatch, sizeOfPatch, u, mappedCell, timeStamp);
    } else {
	interpolateCartesianPatch(offsetOfPatch, sizeOfPatch, u, mappedCell, timeStamp);
    }

    writer->plotPatch(offsetOfPatch, sizeOfPatch, dx, mappedCell, timeStamp);

    delete[] mappedCell;
}


void exahype::plotters::ADERDG2CarpetHDF5::startPlotting(double time) {
	_postProcessing->startPlotting(time);
	writer->startPlotting(time);
}

void exahype::plotters::ADERDG2CarpetHDF5::finishPlotting() {
	_postProcessing->finishPlotting();
	writer->finishPlotting();
}

void exahype::plotters::ADERDG2CarpetHDF5::interpolateCartesianPatch(
  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
  const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
  double *u,
  double *mappedCell,
  double timeStamp
) {
  const int basisSize = writer->basisSize;
  const int solverUnknowns = writer->solverUnknowns;
  const int order = basisSize-1;

  double* interpoland = new double[solverUnknowns];
  //double* mappedCell  = new double[writtenCellIdx.size];
  //double* value       = writtenUnknowns==0 ? nullptr : new double[writtenUnknowns];
  
  assertion(sizeOfPatch(0)==sizeOfPatch(1)); // expressing this is all for squared cells.
  
  dfor(i,basisSize) {
    for (int unknown=0; unknown < solverUnknowns; unknown++) {
      interpoland[unknown] = 0.0;
      dfor(ii,basisSize) { // Gauss-Legendre node indices
        int iGauss = peano::utils::dLinearisedWithoutLookup(ii,order + 1);
        interpoland[unknown] +=
		kernels::equidistantGridProjector1d[order][ii(0)][i(0)] *
		kernels::equidistantGridProjector1d[order][ii(1)][i(1)] *
		#if DIMENSIONS==3
		kernels::equidistantGridProjector1d[order][ii(2)][i(2)] *
		#endif
		u[iGauss * solverUnknowns + unknown];
        assertion3(interpoland[unknown] == interpoland[unknown], offsetOfPatch, sizeOfPatch, iGauss);
      }
    }

    double *value = mappedCell + (DIMENSIONS == 3 ? writer->writtenCellIdx->get(i(2),i(1),i(0),0) : writer->writtenCellIdx->get(i(1),i(0),0));
    //value += writer->writtenCellIdx(i(1),i(0),0); // Transposed position. Correct.

    _postProcessing->mapQuantities(
      offsetOfPatch,
      sizeOfPatch,
      offsetOfPatch + i.convertScalar<double>()* (sizeOfPatch(0)/(order)),
      i,
      interpoland,
      value,
      timeStamp
    );
  }

  delete[] interpoland;
  //if (value!=nullptr)        delete[] value;
}


void exahype::plotters::ADERDG2CarpetHDF5::interpolateSlicedPatch(
  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
  const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
  double *u,
  double *mappedCell,
  double timeStamp
) {
  const int basisSize = writer->basisSize;
  const int solverUnknowns = writer->solverUnknowns;
  const int order = basisSize-1;

  double* interpoland = new double[solverUnknowns];  
  assertion(sizeOfPatch(0)==sizeOfPatch(1)); // expressing this is all for squared cells.
  assertion(writer->slicer!=nullptr); // this is only relevant when the slicer has information
  
  if(writer->slicer->targetDim == 2) {
	// Determine a position ontop the 2d plane
	dvec plane = writer->slicer->project(offsetOfPatch);
	ivec i;
	for(i(1)=0; i(1)<basisSize; i(1)++)
	for(i(0)=0; i(0)<basisSize; i(0)++) {
		dvec pos = plane + writer->slicer->project(i).convertScalar<double>() * (sizeOfPatch(0)/(order));
		
		for (int unknown=0; unknown < solverUnknowns; unknown++) {
			interpoland[unknown] = kernels::interpolate(
				offsetOfPatch.data(),
				sizeOfPatch.data(),
				pos.data(),
				solverUnknowns,
				unknown,
				order,
				u
			);
		}
		
		double *value = mappedCell + writer->writtenCellIdx->get(i(1),i(0),0);
		
		_postProcessing->mapQuantities(
			offsetOfPatch,
			sizeOfPatch,
			pos,
			i,
			interpoland,
			value,
			timeStamp
		);
	}
  } else if(writer->slicer->targetDim == 1) {
	// Determine a position ontop the 1d line
	dvec line = writer->slicer->project(offsetOfPatch);
	ivec i;
	for(i(0)=0; i(0)<basisSize; i(0)++) {
		dvec pos = line + (i.convertScalar<double>())* (sizeOfPatch(0)/(order));
		
		for (int unknown=0; unknown < solverUnknowns; unknown++) {
			interpoland[unknown] = kernels::interpolate(
				offsetOfPatch.data(),
				sizeOfPatch.data(),
				pos.data(),
				solverUnknowns,
				unknown,
				order,
				u
			);
		}
		
		double *value = mappedCell + writer->writtenCellIdx->get(i(0));
		
		_postProcessing->mapQuantities(
			offsetOfPatch,
			sizeOfPatch,
			pos,
			i,
			interpoland,
			value,
			timeStamp
		);
	} 
  } else {
	  throw std::invalid_argument("Unupported target dimension.");
  }
  
  delete[] interpoland;
}



#endif /* HDF5 */
