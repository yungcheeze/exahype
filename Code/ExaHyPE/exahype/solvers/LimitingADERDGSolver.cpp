/*
 * LimitedADERDGSolver.cpp
 *
 *  Created on: 26 Oct 2016
 *      Author: dominic
 */

#include "LimitedADERDGSolver.h"

namespace exahype {
namespace solvers {


tarch::multicore::BooleanSemaphore exahype::solvers::LimitedADERDGSolver::_semaphoreForLimiterHeapAccess;

} /* namespace solvers */
} /* namespace exahype */
