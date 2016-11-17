// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================
#include "exahype/plotters/Plotter.h"
namespace MHDSolver{
  class MHDSolver_Plotter4;

  /**
   * Forward declaration
   */
  class MHDSolver;
}




class MHDSolver::MHDSolver_Plotter4: public exahype::plotters::Plotter::UserOnTheFlyPostProcessing{
  public:
  MHDSolver_Plotter4(MHDSolver&  solver);
  virtual ~MHDSolver_Plotter4();
  void startPlotting(double time) override;
  void finishPlotting() override;
  void mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp) override;
};
