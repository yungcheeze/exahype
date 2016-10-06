/*
 * ADERDGAPosterioriSubcellLimiter.h
 *
 *  Created on: 6 Oct 2016
 *      Author: dominic
 */

#ifndef ADERDGAPOSTERIORISUBCELLLIMITER_H_
#define ADERDGAPOSTERIORISUBCELLLIMITER_H_

#include "CellWiseCoupling.h"

namespace exahype {
namespace solvers {

class ADERDGAPosterioriSubcellLimiter;

} /* namespace solvers */
} /* namespace exahype */

class exahype::solvers::ADERDGAPosterioriSubcellLimiter : public exahype::solvers::CellWiseCoupling {
private:
  bool _isActive;

  const int  _aderdgSolverNumber;
  const int  _finiteVolumesSolverNumber;

public:
  ADERDGAPosterioriSubcellLimiter(int aderdgSolverNumber,int finiteVolumesSolverNumber);
  virtual ~ADERDGAPosterioriSubcellLimiter() {};

  // Disallow copy and assignment
  ADERDGAPosterioriSubcellLimiter(const ADERDGAPosterioriSubcellLimiter& other) = delete;
  ADERDGAPosterioriSubcellLimiter& operator=(const ADERDGAPosterioriSubcellLimiter& other) = delete;

  /**
   * Couple solvers before the solution update of the solvers is performed.
   */
  void coupleSolversBeforeSolutionUpdate(const int cellDescriptionsIndex) override;

  /**
   * Couple solvers after the solution update of the solvers has been performed.
   */
  void coupleSolversAfterSolutionUpdate(const int cellDescriptionsIndex) override;
};

#endif /* ADERDGAPOSTERIORISUBCELLLIMITER_H_ */
