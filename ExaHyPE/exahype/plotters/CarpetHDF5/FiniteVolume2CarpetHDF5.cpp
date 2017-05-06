#include "exahype/plotters/CarpetHDF5/FiniteVolume2CarpetHDF5.h"
#include "exahype/solvers/FiniteVolumesSolver.h"

#include <cstdlib>
#include <stdio.h>
#include <sstream>
#include <memory>


std::string exahype::plotters::FiniteVolume2CarpetHDF5::getIdentifier() {
	return std::string("experimental::CarpetHDF5");
}


// my small C++11 to_string-independent workaround.
template <typename T> std::string toString( T Number ) {
	std::ostringstream ss; ss << Number; return ss.str();
}

typedef tarch::la::Vector<DIMENSIONS, double> dvec;

#ifndef HDF5
/*************************************************************************************************
 * FiniteVolume2CarpetHDF5 Dummy implementation in case HDF5 support is skipped.
 *************************************************************************************************/

exahype::plotters::FiniteVolume2CarpetHDF5::FiniteVolume2CarpetHDF5(
  exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
  const int _ghostLayerWidth) : Device(postProcessing), ghostLayerWidth(_ghostLayerWidth) {
	printf("ERROR: Compile with HDF5, otherwise you cannot use the HDF5 plotter.\n");
	abort();
}

// all other methods are stubs
exahype::plotters::FiniteVolume2CarpetHDF5::~FiniteVolume2CarpetHDF5() {}
void exahype::plotters::FiniteVolume2CarpetHDF5::init(const std::string& filename, int orderPlusOne, int solverUnknowns, int writtenUnknowns, const std::string& select) {}
void exahype::plotters::FiniteVolume2CarpetHDF5::plotPatch(const int cellDescriptionsIndex, const int element) {}
void exahype::plotters::FiniteVolume2CarpetHDF5::plotPatch(const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,double timeStamp) {}
void exahype::plotters::FiniteVolume2CarpetHDF5::startPlotting(double time) {}
void exahype::plotters::FiniteVolume2CarpetHDF5::finishPlotting() {}

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
 * FiniteVolume2CarpetHDF5 non-dummy implementation
 *************************************************************************************************/

exahype::plotters::FiniteVolume2CarpetHDF5::FiniteVolume2CarpetHDF5(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
  const int _ghostLayerWidth) :  Device(postProcessing), ghostLayerWidth(_ghostLayerWidth) { writer = nullptr; }

exahype::plotters::FiniteVolume2CarpetHDF5::~FiniteVolume2CarpetHDF5() {
	if(writer) delete writer;
}

void exahype::plotters::FiniteVolume2CarpetHDF5::init(const std::string& filename, int numberOfCellsPerAxis, int solverUnknowns, int writtenUnknowns, const std::string& select) {
	bool oneFilePerTimestep = true;
	bool allUnknownsInOneFile = true;

	// Determine names of output fields
	char **writtenQuantitiesNames = new char*[writtenUnknowns];
	std::fill_n(writtenQuantitiesNames, writtenUnknowns, nullptr);
	_postProcessing->writtenQuantitiesNames(writtenQuantitiesNames);
	
	const int basisSize = numberOfCellsPerAxis; // apparently this *is* the patchSize
	writer = new exahype::plotters::CarpetHDF5Writer(filename, basisSize, solverUnknowns, writtenUnknowns, select,
		writtenQuantitiesNames, oneFilePerTimestep, allUnknownsInOneFile);	

	// todo at this place: Accept 3D slicing or so, cf. ADERDG2CartesianVTK
	// another nice to have: allow oneFilePerTimestep as a specfile parameter...
}

void exahype::plotters::FiniteVolume2CarpetHDF5::plotPatch(
        const int cellDescriptionsIndex,
        const int element) {
  auto& cellDescription =  exahype::solvers::FiniteVolumesSolver::getCellDescription(cellDescriptionsIndex,element);

  if (cellDescription.getType()==exahype::solvers::FiniteVolumesSolver::CellDescription::Type::Cell) {
    double* solution = DataHeap::getInstance().getData(cellDescription.getSolution()).data();

    plotPatch(
        cellDescription.getOffset(),
        cellDescription.getSize(), solution,
        cellDescription.getTimeStamp());
  }
}

void exahype::plotters::FiniteVolume2CarpetHDF5::plotPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    double* u, /* unknown */
    double timeStamp) {

    double* mappedCell  = new double[writer->writtenCellIdx.size];
    dvec dx = 1./(writer->basisSize) * sizeOfPatch;

    mapCartesianPatch(offsetOfPatch, sizeOfPatch, u, mappedCell, timeStamp);
    writer->plotPatch(offsetOfPatch, sizeOfPatch, dx, mappedCell, timeStamp);

    delete[] mappedCell;
}


void exahype::plotters::FiniteVolume2CarpetHDF5::mapCartesianPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    double* u, /* ingoing unknowns in cell, size numberOfCellsPerAxis^DIMENSIONS */
    double *mappedCell, /* outgoing mapped cell, size writtenUnknowns^DIMENSIONS */
    double timeStamp) {

    const int basisSize = writer->basisSize; // == numberOfCellsPerAxis == nVar*patchSize^DIMENSIONS
    const int writtenUnknowns = writer->writtenUnknowns;
    const int solverUnknowns = writer->solverUnknowns;

    dfor(i,basisSize) {
        double *value = mappedCell;
        #if DIMENSIONS==3
        value += writer->writtenCellIdx(i(2),i(1),i(0),0);
        #else
        value += writer->writtenCellIdx(i(1),i(0),0);
        #endif

	_postProcessing->mapQuantities(
          offsetOfPatch,
          sizeOfPatch,
          offsetOfPatch + i.convertScalar<double>()*(sizeOfPatch(0)/basisSize),
          i,
          u + peano::utils::dLinearisedWithoutLookup(i+ghostLayerWidth, basisSize+2*ghostLayerWidth) * solverUnknowns,
          value,
          timeStamp
        );
    }
}


void exahype::plotters::FiniteVolume2CarpetHDF5::startPlotting(double time) {
	_postProcessing->startPlotting(time);
	writer->startPlotting(time);
}

void exahype::plotters::FiniteVolume2CarpetHDF5::finishPlotting() {
	_postProcessing->finishPlotting();
	writer->finishPlotting();
}

#endif /* HDF5 */
