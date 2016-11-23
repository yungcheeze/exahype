/**
 * This is the fancy FORTRAN and C calling kernel generated stuff.
 *
 **/
#include "MHDSolver.h"

// includes fortran and c kernels
#include "kernels/aderdg/generic/Kernels.h"
#include "kernels/KernelUtils.h"

// std::abs
#include <climits>
#include <cstdlib>


MHDSolver::MHDSolver::MHDSolver(double maximumMeshSize,exahype::solvers::Solver::TimeStepping timeStepping,std::vector<std::string>& cmdlineargs, exahype::Parser::ParserView constants):
  exahype::solvers::ADERDGSolver("MHDSolver",nVar /* numberOfUnknowns */,nParams /* numberOfParameters */,order + 1  /* nodesPerCoordinateAxis */,maximumMeshSize,timeStepping) {
  init(cmdlineargs, constants);
}


int globalspaceTimePredictorCounter(0);

void MHDSolver::MHDSolver::spaceTimePredictor(double* lQhbnd,double* lFhbnd,double** tempSpaceTimeUnknowns,double** tempSpaceTimeFluxUnknowns,double* tempUnknowns,double* tempFluxUnknowns,double* tempStateSizedVectors,const double* const luh,const tarch::la::Vector<DIMENSIONS,double>& dx,const double dt) {
  const int basisSize = order + 1;

  globalspaceTimePredictorCounter++;
  // Fortran comparision: Output arrays
  kernels::idx4 idx_lQhbnd(2 * DIMENSIONS, basisSize, basisSize, nVar);
  kernels::idx4 idx_lFhbnd(2 * DIMENSIONS, basisSize, basisSize, nVar);
  const int size = 2 * DIMENSIONS*basisSize* basisSize*nVar;
  double fort_lQhbnd[size];
  double fort_lFhbnd[size];

  kernels::aderdg::generic::fortran::spaceTimePredictorNonlinear<MHDSolver>(*this,fort_lQhbnd,fort_lFhbnd,tempSpaceTimeUnknowns,tempSpaceTimeFluxUnknowns,tempUnknowns,tempFluxUnknowns,tempStateSizedVectors,luh,dx,dt);
  kernels::aderdg::generic::c::spaceTimePredictorNonlinear<MHDSolver>(*this,lQhbnd,lFhbnd,tempSpaceTimeUnknowns,tempSpaceTimeFluxUnknowns,tempUnknowns,tempFluxUnknowns,tempStateSizedVectors,luh,dx,dt);

  const double accepterror = 1e-10;
  int errors = 0;
  double minError = 0, maxError = 1234;
  for(int i0=0; i0<2*DIMENSIONS; i0++) {
  for(int i1=0; i1<basisSize; i1++) {
  for(int i2=0; i2<basisSize; i2++) {
  for(int i3=0; i3<nVar; i3++) {
	int i = idx_lQhbnd(i0,i1,i2,i3);
	double errorQh = std::abs(fort_lQhbnd[i] - lQhbnd[i]);
	double errorFh = std::abs(fort_lFhbnd[i] - lFhbnd[i]);
	if(errorQh > accepterror || errorFh > accepterror) {
		//printf("%d. idx(%i,%i,%i,%i) lQCPP[%i] = %+.20e lQhFORT[%i] = %+.20e error=%+e\n", globalspaceTimePredictorCounter, i0,i1,i2,i3, i, lQhbnd[i], i, fort_lQhbnd[i], errorQh);
		//printf("%d. idx(%i,%i,%i,%i) lFCPP[%i] = %+.20e lFhFORT[%i] = %+.20e error=%+e\n", globalspaceTimePredictorCounter, i0,i1,i2,i3, i, lFhbnd[i], i, fort_lFhbnd[i], errorFh);
		errors++;
	}
	// log for statistics
	minError = std::min(minError, errorQh);
	minError = std::min(minError, errorFh);
	maxError = std::max(maxError, errorQh);
	maxError = std::max(maxError, errorFh);
  }
  }
  }
  }
  if(errors>0) {
    printf("spaceTimePredictor: %d/%d errors > %e in spaceTimePredictor (minError=%e, maxError=%e)\n", errors, size, accepterror, minError, maxError);
    exit(-1);
  }
}



void MHDSolver::MHDSolver::solutionUpdate(double* luh,const double* const lduh,const double dt) {
  kernels::aderdg::generic::c::solutionUpdate(luh,lduh,dt,getNumberOfVariables(),getNumberOfParameters(),getNodesPerCoordinateAxis());
}



void MHDSolver::MHDSolver::volumeIntegral(double* lduh,const double* const lFhi,const tarch::la::Vector<DIMENSIONS,double>& dx) {
  kernels::aderdg::generic::c::volumeIntegralNonlinear(lduh,lFhi,dx,getNumberOfVariables(),getNumberOfParameters(),getNodesPerCoordinateAxis());
}



void MHDSolver::MHDSolver::surfaceIntegral(double* lduh,const double* const lFhbnd,const tarch::la::Vector<DIMENSIONS,double>& dx) {
  kernels::aderdg::generic::c::surfaceIntegralNonlinear(lduh,lFhbnd,dx,getNumberOfVariables(),getNodesPerCoordinateAxis());
}

int globalRiemannSolverCounter(0);

void MHDSolver::MHDSolver::riemannSolver(double* FL,double* FR,const double* const QL,const double* const QR,double* tempFaceUnknownsArray,double** tempStateSizedVectors,double** tempStateSizedSquareMatrices,const double dt,const int normalNonZeroIndex) {
  assertion2(normalNonZeroIndex>=0,dt,normalNonZeroIndex);
  assertion2(normalNonZeroIndex<DIMENSIONS,dt,normalNonZeroIndex);
  const int basissize = order + 1;
  globalRiemannSolverCounter++;
#define FORTRANCOMPARISON
#ifdef FORTRANCOMPARISON
  // First call Fortran as we know it works and doesn't corrupt data.
  // output quantities
  double fort_FL[basissize*basissize*nVar];
  double fort_FR[basissize*basissize*nVar];
  std::memcpy(fort_FL, FL, basissize*basissize*nVar * sizeof(double));
  std::memcpy(fort_FR, FR, basissize*basissize*nVar * sizeof(double));
  // input quantities. May be copied for security. Not used do far.
  double fort_QL[basissize*basissize*nVar];
  double fort_QR[basissize*basissize*nVar];
  // temporary storage for computation
  double *fort_tempFaceUnknownsArray = NULL; // This is never used.
  const int numTempStateSizedVectors = 5; // cf. ADERDGSolver.h, line 738ff.
  double *fort_tempStateSizedVectors[numTempStateSizedVectors];
  double v1[numTempStateSizedVectors][nVar];
  for(int i=0; i<numTempStateSizedVectors; i++) {
    fort_tempStateSizedVectors[i] = v1[i];
  }
  double **fort_tempStateSizedSquareMatrices = NULL; // Never used

  kernels::aderdg::generic::fortran::riemannSolverNonlinear<MHDSolver>(*this,fort_FL,fort_FR,QL,QR,tempFaceUnknownsArray,tempStateSizedVectors,tempStateSizedSquareMatrices,dt,normalNonZeroIndex);
#endif
  // then call C on the actual data
  kernels::aderdg::generic::c::riemannSolverNonlinear<MHDSolver>(*this,FL,FR,QL,QR,tempFaceUnknownsArray,tempStateSizedVectors,tempStateSizedSquareMatrices,dt,normalNonZeroIndex);
#ifdef FORTRANCOMPARISON
  // now compare data
  const double accepterror = 1e-10;
  bool error = false;
  kernels::idx3 idx_FQLR(basissize, basissize, nVar);
  for (int n = 0; n < basissize; n++) {
    for (int j = 0; j < basissize; j++) {
        for (int k = 0; k < nVar - 0; k++) {
		int i = idx_FQLR(n, j, k);
		double errorFL = std::abs(fort_FL[i] - FL[i]);
		double errorFR = std::abs(fort_FR[i] - FR[i]);
		if(errorFL > accepterror || errorFR > accepterror) {
			printf("idx(%d,%d,%d) = [%d]\n", n,j,k,i);
			printf("%d. FL[%d]=%+e fort_FL[%d]=%+e errorFL=%e\n", globalRiemannSolverCounter, i, FL[i], i, fort_FL[i], errorFL);
			printf("%d. FR[%d]=%+e fort_FR[%d]=%+e errorFR=%e\n", globalRiemannSolverCounter, i, FR[i], i, fort_FR[i], errorFR);
			error = true;
		}
		
	}
    }
  }
  if(!error) {
   // printf("%d. correct\n", globalRiemannSolverCounter);
  } else {
    printf("Stopping since error\n");
    exit(-1);
  }
#endif
}



void MHDSolver::MHDSolver::boundaryConditions(double* fluxOut,double* stateOut,const double* const fluxIn,const double* const stateIn,const tarch::la::Vector<DIMENSIONS,double>& cellCentre,const tarch::la::Vector<DIMENSIONS,double>& cellSize,const double t,const double dt,const int faceIndex,const int normalNonZero) {
  kernels::aderdg::generic::c::boundaryConditions<MHDSolver>(*this,fluxOut,stateOut,fluxIn,stateIn,cellCentre,cellSize,t,dt,faceIndex,normalNonZero);
}



double MHDSolver::MHDSolver::stableTimeStepSize(const double* const luh,double* tempEigenvalues,const tarch::la::Vector<DIMENSIONS,double>& dx) {
  double d = kernels::aderdg::generic::c::stableTimeStepSize<MHDSolver>(*this,luh,tempEigenvalues,dx);
  return d;
}



void MHDSolver::MHDSolver::solutionAdjustment(double *luh,const tarch::la::Vector<DIMENSIONS,double>& center,const tarch::la::Vector<DIMENSIONS,double>& dx,double t,double dt) {
  kernels::aderdg::generic::c::solutionAdjustment<MHDSolver>(*this,luh,center,dx,t,dt);
}



void MHDSolver::MHDSolver::faceUnknownsProlongation(double* lQhbndFine,double* lFhbndFine,const double* lQhbndCoarse,const double* lFhbndCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1,int>& subfaceIndex) {
  kernels::aderdg::generic::c::faceUnknownsProlongation(lQhbndFine,lFhbndFine,lQhbndCoarse,lFhbndCoarse,coarseGridLevel,fineGridLevel,subfaceIndex,getNumberOfVariables(),getNodesPerCoordinateAxis());
}



void MHDSolver::MHDSolver::faceUnknownsRestriction(double* lQhbndCoarse,double* lFhbndCoarse,const double* lQhbndFine,const double* lFhbndFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS-1,int>& subfaceIndex) {
  kernels::aderdg::generic::c::faceUnknownsRestriction(lQhbndCoarse,lFhbndCoarse,lQhbndFine,lFhbndFine,coarseGridLevel,fineGridLevel,subfaceIndex,getNumberOfVariables(),getNodesPerCoordinateAxis());
}



void MHDSolver::MHDSolver::volumeUnknownsProlongation(double* luhFine,const double* luhCoarse,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) {
  kernels::aderdg::generic::c::volumeUnknownsProlongation(luhFine,luhCoarse,coarseGridLevel,fineGridLevel,subcellIndex,getNumberOfVariables(),getNodesPerCoordinateAxis());
}



void MHDSolver::MHDSolver::volumeUnknownsRestriction(double* luhCoarse,const double* luhFine,const int coarseGridLevel,const int fineGridLevel,const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) {
  kernels::aderdg::generic::c::volumeUnknownsRestriction(luhCoarse,luhFine,coarseGridLevel,fineGridLevel,subcellIndex,getNumberOfVariables(),getNodesPerCoordinateAxis());
}