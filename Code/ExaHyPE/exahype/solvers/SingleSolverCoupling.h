/*
 * SingleSolverCoupling.h
 *
 *  Created on: 6 Oct 2016
 *      Author: dominic
 */

#ifndef SingleSolverCoupling_H_
#define SingleSolverCoupling_H_

#include "CellWiseCoupling.h"

namespace exahype {
namespace solvers {

class SingleSolverCoupling;

} /* namespace solvers */
} /* namespace exahype */

/**
 * This class simply performs a solution update for the wrapped solver.
 */
class exahype::solvers::SingleSolverCoupling : public exahype::solvers::CellWiseCoupling {
private:
  /**
   * Element index solver exahype::solvers::RegisteredSolvers.
   */
  const int _solverNumber;
public:
  SingleSolverCoupling(int solverNumber);
  virtual ~SingleSolverCoupling() {};

  // Disallow copy and assignment
  SingleSolverCoupling(const SingleSolverCoupling& other) = delete;
  SingleSolverCoupling& operator=(const SingleSolverCoupling& other) = delete;

  /**
   * Simply calls the method setInitialConditions() on the wrapped solver.
   */
  void coupleFirstTime(
      const int cellDescriptionsIndex,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) override;

  /**
   * Simply calls the method updateSolution() on the wrapped solver.
   */
  void couple(
      const int cellDescriptionsIndex,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) override;
};

#endif /* SingleSolverCoupling_H_ */
