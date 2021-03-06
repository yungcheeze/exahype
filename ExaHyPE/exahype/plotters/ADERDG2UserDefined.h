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
 
#ifndef _EXAHYPE_PLOTTERS_ADERDG_2_USER_DEFINED_H_
#define _EXAHYPE_PLOTTERS_ADERDG_2_USER_DEFINED_H_

#include "exahype/plotters/Plotter.h"

namespace exahype {
  namespace plotters {
    class ADERDG2UserDefined;
}
}

/**
 * Device for realising user defined plotters.
 */
class exahype::plotters::ADERDG2UserDefined: public exahype::plotters::Plotter::Device {
 protected:
  std::string   _filename;
  int           _order;
  int           _variables;
  int           _writtenVariables;
  std::string   _select;

 public:
  static std::string getIdentifier();

  ADERDG2UserDefined();
  virtual ~ADERDG2UserDefined();

  void init(const std::string& filename, int orderPlusOne, int solverUnknowns, int writtenUnknowns, const std::string& select) override;

  void plotPatch(const int cellDescriptionsIndex, const int element) override;

  virtual void plotPatch(
      const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
      const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
      double timeStamp) = 0;

  virtual void startPlotting( double time) = 0;
  virtual void finishPlotting() = 0;
};


#endif
