#include "exahype/tests/kernels/fortran/GenericEulerKernelTest.h"

#include "kernels/aderdg/generic/Kernels.h"

using std::cout;
using std::endl;

#if DIMENSIONS == 3

void exahype::tests::fortran::GenericEulerKernelTest::testPDEFluxes() {
  cout << "Test PDE-related functions, DIM=3" << endl;

}  // testPDEFluxes


void exahype::tests::fortran::GenericEulerKernelTest::testSolutionUpdate() {
  cout << "Test solution update, ORDER=3, DIM=3" << endl;

}  // testSolutionUpdate


void exahype::tests::fortran::GenericEulerKernelTest::testSurfaceIntegral() {
  cout << "Test surface integral, ORDER=3, DIM=3" << endl;

}  // testSurfaceIntegral

void exahype::tests::fortran::GenericEulerKernelTest::testRiemannSolver() {

  cout << "Test Riemann Solver (Rusanov), ORDER=3, DIM=3" << endl;
}  // testRiemannSolver

void exahype::tests::fortran::GenericEulerKernelTest::testVolumeIntegral() {
  cout << "Test volume integral, ORDER=3, DIM=3" << endl;

}  // testVolumeIntegral

void exahype::tests::fortran::GenericEulerKernelTest::testSpaceTimePredictor() {

  cout << "Test space time predictor, ORDER=3, DIM=3" << endl;
}  // testSpaceTimePredictor

#endif  // DIMENSIONS==3
