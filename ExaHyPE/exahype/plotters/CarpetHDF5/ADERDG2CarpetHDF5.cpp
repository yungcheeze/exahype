#include "exahype/plotters/CarpetHDF5/ADERDG2CarpetHDF5.h"

#include <cstdlib>
#include <stdio.h>
#include <sstream>
#include <memory>


std::string exahype::plotters::ADERDG2CarpetHDF5::getIdentifier() {
	return std::string("experimental::CarpetHDF5");
}


// my small C++11 to_string-independent workaround.
template <typename T> std::string toString( T Number ) {
	std::ostringstream ss; ss << Number; return ss.str();
}

typedef tarch::la::Vector<DIMENSIONS, double> dvec;

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
	bool allUnknownsInOneFile = false;

	// Determine names of output fields
	char **writtenQuantitiesNames = new char*[writtenUnknowns];
	std::fill_n(writtenQuantitiesNames, writtenUnknowns, nullptr);
	_postProcessing->writtenQuantitiesNames(writtenQuantitiesNames);

	writer = new exahype::plotters::CarpetHDF5Writer(filename, basisSize, solverUnknowns, writtenUnknowns, select,
		writtenQuantitiesNames, oneFilePerTimestep, allUnknownsInOneFile);	

	// todo at this place: Accept 3D slicing or so, cf. ADERDG2CartesianVTK
	// another nice to have: allow oneFilePerTimestep as a specfile parameter...
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
    const int basisSize = writer->basisSize;
    const int order = basisSize - 1;


    // TODO: if we knew that plotting would be serial, we could move *mappedCell to a class property.
    double* mappedCell  = new double[writer->writtenCellIdx.size];

    dvec dx = 1./order * sizeOfPatch;

    interpolateCartesianPatch(offsetOfPatch, sizeOfPatch, u, mappedCell, timeStamp);
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

    double *value = mappedCell;
    #if DIMENSIONS==3
    value += writer->writtenCellIdx(i(2),i(1),i(0),0);
    #else
    value += writer->writtenCellIdx(i(1),i(0),0);
    #endif
    assertion(sizeOfPatch(0)==sizeOfPatch(1));
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




#endif /* HDF5 */
