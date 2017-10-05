/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/

#include "exahype/solvers/TemporaryVariables.h"

#include "exahype/solvers/Solver.h"
#include "exahype/solvers/LimitingADERDGSolver.h"

#include <algorithm>
#include <mm_malloc.h> //g++
#include <cstring> //memset

#include "tarch/Assertions.h"

exahype::solvers::PredictionTemporaryVariables::PredictionTemporaryVariables() {}

double* exahype::solvers::allocateArray(int size, bool align) {
  //don't need to allocate anything
  if(size < 1) {
    return nullptr;
  }
  double* array;
  #ifdef ALIGNMENT
  if(align) {
    array = ((double *) _mm_malloc(sizeof(double)*size, ALIGNMENT));
  } else {  
    array = new double[size]();
  }
  #else
    array = new double[size]();
  #endif
  std::memset(array, 0, sizeof(double)*size);
  
  return array;
}

void exahype::solvers::freeArray(double*& array, bool align) {
  //don't need to delete
  if(array == nullptr) {
     return;
  }

  #ifdef ALIGNMENT
  if(align) {
  _mm_free(array);
  } else {
    delete[] array;
  }
  #else
  delete[] array;
  #endif
  array = nullptr;
}

void exahype::solvers::initialiseTemporaryVariables(exahype::solvers::PredictionTemporaryVariables& temporaryVariables) {
  assertion(temporaryVariables._tempSpaceTimeUnknowns    ==nullptr);
  assertion(temporaryVariables._tempSpaceTimeFluxUnknowns==nullptr);
  assertion(temporaryVariables._tempUnknowns             ==nullptr);
  assertion(temporaryVariables._tempFluxUnknowns         ==nullptr);
  assertion(temporaryVariables._tempPointForceSources    ==nullptr);

  int numberOfSolvers        = exahype::solvers::RegisteredSolvers.size();
  temporaryVariables._tempSpaceTimeUnknowns     = new double**[numberOfSolvers]; // == lQi, rhs
  temporaryVariables._tempSpaceTimeFluxUnknowns = new double**[numberOfSolvers]; // == lFi, gradQ
  temporaryVariables._tempUnknowns              = new double* [numberOfSolvers]; // == lQhi
  temporaryVariables._tempFluxUnknowns          = new double* [numberOfSolvers]; // == lFhi
  temporaryVariables._tempPointForceSources     = new double* [numberOfSolvers];

  exahype::solvers::ADERDGSolver* aderdgSolver = nullptr;

  int solverNumber=0;
  for (auto solver : exahype::solvers::RegisteredSolvers) {
    switch( solver->getType() ) {
    case exahype::solvers::Solver::Type::ADERDG:
      aderdgSolver = static_cast<exahype::solvers::ADERDGSolver*>(solver);
      break;
    case exahype::solvers::Solver::Type::LimitingADERDG:
      aderdgSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver().get();
      break;
    default:
      aderdgSolver = nullptr;
      break;
    }

    if (aderdgSolver!=nullptr) {
      const bool align = aderdgSolver->alignTempArray();
      
      temporaryVariables._tempSpaceTimeUnknowns[solverNumber] = new double*[2];
      for (int i=0; i<2; ++i) { // max; see spaceTimePredictorNonlinear
        temporaryVariables._tempSpaceTimeUnknowns[solverNumber][i] = allocateArray(aderdgSolver->getTempSpaceTimeUnknownsSize(), align);
      }
      //
      temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber] = new double*[2];
      for (int i=0; i<2; ++i) { // max; see spaceTimePredictorNonlinear
        temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber][i] = allocateArray(aderdgSolver->getTempSpaceTimeFluxUnknownsSize(), align);
      }
      //
      temporaryVariables._tempUnknowns    [solverNumber]      = allocateArray(aderdgSolver->getTempUnknownsSize(), align);
      //
      temporaryVariables._tempFluxUnknowns[solverNumber]      = allocateArray(aderdgSolver->getTempFluxUnknownsSize(), align);
      //
      temporaryVariables._tempPointForceSources    [solverNumber] = allocateArray(aderdgSolver->getTempPointForceSourcesSize(), align);
    } else {
      temporaryVariables._tempSpaceTimeUnknowns    [solverNumber] = nullptr;
      temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber] = nullptr;
      temporaryVariables._tempUnknowns             [solverNumber] = nullptr;
      temporaryVariables._tempFluxUnknowns         [solverNumber] = nullptr;
      temporaryVariables._tempPointForceSources    [solverNumber] = nullptr;
    }

    ++solverNumber;
  }
}

void exahype::solvers::deleteTemporaryVariables(exahype::solvers::PredictionTemporaryVariables& temporaryVariables) {
  if (temporaryVariables._tempSpaceTimeUnknowns!=nullptr) {
    assertion(temporaryVariables._tempSpaceTimeFluxUnknowns!=nullptr);
    //assertion(temporaryVariables._tempUnknowns             !=nullptr);
    //assertion(temporaryVariables._tempFluxUnknowns         !=nullptr);
    assertion(temporaryVariables._tempPointForceSources    !=nullptr);

    int solverNumber=0;
    exahype::solvers::ADERDGSolver* aderdgSolver = nullptr;

    for (auto solver : exahype::solvers::RegisteredSolvers) {
      switch( solver->getType() ) {
      case exahype::solvers::Solver::Type::ADERDG:
        aderdgSolver = static_cast<exahype::solvers::ADERDGSolver*>(solver);
        break;
      case exahype::solvers::Solver::Type::LimitingADERDG:
        aderdgSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver().get();
        break;
      default:
        aderdgSolver = nullptr;
        break;
      }

      if (aderdgSolver!=nullptr) {
          const bool align = aderdgSolver->alignTempArray();
          //
          for (int i=0; i<2; ++i) {
            freeArray(temporaryVariables._tempSpaceTimeUnknowns[solverNumber][i], align);
          }
          delete[] temporaryVariables._tempSpaceTimeUnknowns[solverNumber];
          temporaryVariables._tempSpaceTimeUnknowns[solverNumber] = nullptr;
          //
          for (int i=0; i<2; ++i) {
            freeArray(temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber][i], align);
          }
          delete[] temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber];
          temporaryVariables._tempSpaceTimeFluxUnknowns[solverNumber] = nullptr;
          //
          freeArray(temporaryVariables._tempUnknowns[solverNumber], align);
          //
          freeArray(temporaryVariables._tempFluxUnknowns[solverNumber], align);
          //
          freeArray(temporaryVariables._tempPointForceSources[solverNumber], align);
      }

      ++solverNumber;
    }

    delete[] temporaryVariables._tempSpaceTimeUnknowns;
    delete[] temporaryVariables._tempSpaceTimeFluxUnknowns;
    delete[] temporaryVariables._tempUnknowns;
    delete[] temporaryVariables._tempFluxUnknowns;
    delete[] temporaryVariables._tempPointForceSources;
    temporaryVariables._tempSpaceTimeUnknowns     = nullptr;
    temporaryVariables._tempSpaceTimeFluxUnknowns = nullptr;
    temporaryVariables._tempUnknowns              = nullptr;
    temporaryVariables._tempFluxUnknowns          = nullptr;
    temporaryVariables._tempPointForceSources     = nullptr;
  }
}

exahype::solvers::MergingTemporaryVariables::MergingTemporaryVariables() {}

void exahype::solvers::initialiseTemporaryVariables(exahype::solvers::MergingTemporaryVariables& temporaryVariables) {
  assertion(temporaryVariables._tempFaceUnknowns            ==nullptr);

  int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
  temporaryVariables._tempFaceUnknowns             = new double**[numberOfSolvers];

  int solverNumber=0;
  for (auto solver : exahype::solvers::RegisteredSolvers) {
    int numberOfFaceUnknowns = 0;
    int lengthOfFaceUnknowns = 0;
    bool alignTempArray      = false;

    switch (solver->getType()) {
      case exahype::solvers::Solver::Type::ADERDG:
        numberOfFaceUnknowns       = 3; // See exahype::solvers::ADERDGSolver::applyBoundaryConditions
        lengthOfFaceUnknowns       =
            static_cast<exahype::solvers::ADERDGSolver*>(solver)->getBndFaceSize(); // == getDataPerFace() + eventual padding
        alignTempArray = static_cast<exahype::solvers::ADERDGSolver*>(solver)->alignTempArray();
        break;
      case exahype::solvers::Solver::Type::LimitingADERDG:
        // Needs the same temporary variables as the normal ADER-DG scheme.
        numberOfFaceUnknowns       = 3;
        lengthOfFaceUnknowns       = std::max(
            static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver()->getBndFaceSize(), // == getDataPerFace() + eventual padding
            static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiter()->getDataPerPatchFace() );
        alignTempArray = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver()->alignTempArray();
        break;
      case exahype::solvers::Solver::Type::FiniteVolumes:
        numberOfFaceUnknowns = 2; // See exahype::solvers::FiniteVolumesSolver::mergeWithBoundaryData
        lengthOfFaceUnknowns =
            static_cast<exahype::solvers::FiniteVolumesSolver*>(solver)->getDataPerPatchFace();
        break;
    }
    //
    temporaryVariables._tempFaceUnknowns[solverNumber] = nullptr;
    if (numberOfFaceUnknowns>0) {
      temporaryVariables._tempFaceUnknowns[solverNumber] = new double*[numberOfFaceUnknowns];
      for (int i=0; i<numberOfFaceUnknowns; ++i) { // see ADERDGSolver::applyBoundaryConditions(...)
        #ifdef ALIGNMENT
        if (alignTempArray) {
          temporaryVariables._tempFaceUnknowns[solverNumber][i] =
              (double*) _mm_malloc( sizeof(double)*lengthOfFaceUnknowns, ALIGNMENT );
        } else {
          temporaryVariables._tempFaceUnknowns[solverNumber][i] = new double[lengthOfFaceUnknowns](); //initialized to 0 to ensure padding is initialized if existing
        }
        #else
        temporaryVariables._tempFaceUnknowns[solverNumber][i] = new double[lengthOfFaceUnknowns](); //initialized to 0 to ensure padding is initialized if existing
        #endif
      }
    }

    ++solverNumber;
  }
}

void exahype::solvers::deleteTemporaryVariables(exahype::solvers::MergingTemporaryVariables& temporaryVariables) {
  if (temporaryVariables._tempFaceUnknowns!=nullptr) {
    int solverNumber=0;
    for (auto solver : exahype::solvers::RegisteredSolvers) {
      int numberOfFaceUnknowns = 0;
      bool alignTempArray      = false;

      switch (solver->getType()) {
        case exahype::solvers::Solver::Type::ADERDG:
          numberOfFaceUnknowns = 3; // See exahype::solvers::ADERDGSolver::applyBoundaryConditions
          alignTempArray       = static_cast<exahype::solvers::ADERDGSolver*>(solver)->alignTempArray();
          break;
        case exahype::solvers::Solver::Type::LimitingADERDG:
          numberOfFaceUnknowns = 3; // See exahype::solvers::ADERDGSolver::applyBoundaryConditions
          alignTempArray       = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver()->alignTempArray();
          break;
        case exahype::solvers::Solver::Type::FiniteVolumes:
          numberOfFaceUnknowns = 2; // See exahype::solvers::FiniteVolumesSolver::mergeWithBoundaryData
          break;
      }

      if (numberOfFaceUnknowns>0) {
        for (int i=0; i<numberOfFaceUnknowns; ++i) { // see ADERDGSolver::applyBoundaryConditions(...)
          #ifdef ALIGNMENT
          if (alignTempArray) {
            _mm_free( temporaryVariables._tempFaceUnknowns[solverNumber][i] );
          } else {
            delete[] temporaryVariables._tempFaceUnknowns[solverNumber][i];
          }
          #else
          delete[] temporaryVariables._tempFaceUnknowns[solverNumber][i];
          #endif
        }
        delete[] temporaryVariables._tempFaceUnknowns[solverNumber];
        temporaryVariables._tempFaceUnknowns[solverNumber] = nullptr;
      }

      ++solverNumber;
    }

    delete[] temporaryVariables._tempFaceUnknowns;
    temporaryVariables._tempFaceUnknowns             = nullptr;
  }
}
