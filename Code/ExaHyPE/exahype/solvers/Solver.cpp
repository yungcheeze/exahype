#include "exahype/solvers/Solver.h"


std::vector<exahype::solvers::Solver*>  exahype::solvers::RegisteredSolvers;


exahype::solvers::Solver::Solver(const std::string& identifier, Type type, int kernelNumber, int numberOfVariables, int nodesPerCoordinateAxis):
          _identifier(identifier),
          _type(type),
          _kernelNumber(kernelNumber),
          _numberOfVariables(numberOfVariables),
          _nodesPerCoordinateAxis(nodesPerCoordinateAxis) {
}


exahype::solvers::Solver::Type exahype::solvers::Solver::getType() const {
  return _type;
}


int exahype::solvers::Solver::getNumberOfVariables() const {
  return _numberOfVariables;
}


int exahype::solvers::Solver::getNodesPerCoordinateAxis() const {
  return _nodesPerCoordinateAxis;
}


std::string exahype::solvers::Solver::getIdentifier() const {
  return _identifier;
}



// @todo Dominic Charrier
// Just to have the definition headers ready for the sub classes.
// I will delete them later on.

//void exahype::solvers::Solver::solutionUpdate(
//    double * luh,
//    const double * const lduh,
//    const tarch::la::Vector<DIMENSIONS,double> dx,
//    const double dt
//) {
//     const double size  [2] = { fineGridVerticesEnumerator.getCellSize()  [0], fineGridVerticesEnumerator.getCellSize()  [1]};
//aderdg::solutionUpdate<DIMENSIONS>(
//    luhOld,
//    lduh,
//    size,
//    _localState.getMaxTimeStepSize(),
//    nvar,
//    basisSize);
//}
//
//void exahype::solvers::Solver::volumeIntegral(
//    double * lduh,
//    const double * const lFhi,
//    const tarch::la::Vector<DIMENSIONS,double> dx
//) {
// aderdg::volumeIntegral<DIMENSIONS>(
//               lduh,
//               lFhi,
//              size);
//}
//
//void exahype::solvers::Solver::surfaceIntegral(
//    double * lduh,
//    const double * const lFhbnd,
//    const tarch::la::Vector<DIMENSIONS,double> dx
//) {
//
//    @todo Dominic Etienne Charrier
//    Not sure where to place these lines.
//    const int dofStartIndexLeft  = EXAHYPE_FACE_LEFT  * numberOfFaceDof;
//    const int dofStartIndexRight = EXAHYPE_FACE_RIGHT * numberOfFaceDof;
//    const int dofStartIndexFront = EXAHYPE_FACE_FRONT * numberOfFaceDof;
//    const int dofStartIndexBack  = EXAHYPE_FACE_BACK  * numberOfFaceDof;
//
//    double * lFhLeft  = &(DataHeap::getInstance().getData(p->getFluctuation())[dofStartIndexLeft ]._persistentRecords._u);
//    double * lFhRight = &(DataHeap::getInstance().getData(p->getFluctuation())[dofStartIndexRight]._persistentRecords._u);
//    double * lFhFront = &(DataHeap::getInstance().getData(p->getFluctuation())[dofStartIndexFront]._persistentRecords._u);
//    double * lFhBack  = &(DataHeap::getInstance().getData(p->getFluctuation())[dofStartIndexBack ]._persistentRecords._u);
//    aderdg::surfaceIntegral(
//        lduh,
//        size,
//        nvar,
//        basisSize,
//        lFhLeft,
//        lFhRight,
//        lFhFront,
//        lFhBack);
//}
//
//void exahype::solvers::Solver::riemannSolver(
//    double * FL,
//    double * FR,
//    const double * const QL,
//    const double * const QR,
//    const double dt,
//    const double hFace,
//    const tarch::la::Vector<DIMENSIONS,double> n
//) {
//
//}
//
//void exahype::solvers::Solver::spaceTimePredictor(
//    double * lQi,
//    double * lFi,
//    double * lQhi,
//    double * lFhi,
//    double * lQhbnd,
//    double * lFhbnd,
//    const double * const luh,
//    const tarch::la::Vector<DIMENSIONS,double> dx,
//    const double dt
//) {
//  // todo temp variables
//}
//
//void exahype::solvers::Solver::initialValues(
//    double * luh,
//    const tarch::la::Vector<DIMENSIONS,double> center,
//    const tarch::la::Vector<DIMENSIONS,double> dx,
//) {
//
//}
//
//void exahype::solvers::Solver::plot(
//    const double * luh,
//    const tarch::la::Vector<DIMENSIONS,double> center,
//    const tarch::la::Vector<DIMENSIONS,double> dx
//) {
//   exahype::aderdg::io::exportToVTK<DIMENSIONS>(_vtkWriter,_vertexWriter,_cellWriter,_vertexValueWriter,luh,center,size,nvar,basisSize);
//}
