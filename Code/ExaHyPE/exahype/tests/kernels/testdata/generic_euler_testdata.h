#ifndef _EXAHYPE_TESTS_TESTDATA_GENERIC_EULER_TESTDATA_H_
#define _EXAHYPE_TESTS_TESTDATA_GENERIC_EULER_TESTDATA_H_

namespace exahype {
namespace tests {
namespace testdata {
namespace generic_euler {

#ifdef Dim2

namespace testPDEFluxes {
extern const double f[5];
extern const double g[5];
}  // namespace testPDEFluxes

namespace testSolutionUpdate {
extern const double lduh[80];
extern const double luh[80];
}  // namespace testSolutionUpdate

namespace testSurfaceIntegral {
extern const double lduh_in[80];
extern const double lduh_out[80];
}  // namespace testSurfaceIntegral

namespace testVolumeIntegralNonlinear {
extern const double lFhi[160];
extern const double lduh_1[80];
extern const double lduh_2[80];
}  // namespace testVolumeIntegralNonlinear

namespace testVolumeIntegralNonlinear {
extern const double lduh_1[80];
extern const double lduh_2[80];
}  // namespace testVolumeIntegralNonlinear

namespace testSpaceTimePredictorNonlinear {
extern const double lQhi[80];
extern const double lFhi[160];
extern const double lQhbnd[80];
extern const double lFhbnd[80];
}  // namespace testSpaceTimePredictorNonlinear

#endif  // Dim2

#ifdef Dim3

namespace testPDEFluxes {
extern const double f[5];
extern const double g[5];
extern const double h[5];
}  // namespace testPDEFluxes

namespace testVolumeIntegralLinear {
extern const double lduh[320];
}  // namespace testVolumeIntegralLinear

namespace testVolumeIntegralNonlinear {
extern const double lduh[320];
}  // namespace testVolumeIntegralNonlinear

namespace testSurfaceIntegral {
extern const double lduh[320];
}  // namesapce testSurfaceIntegral

namespace testRiemannSolver {
extern const double QL[80];
extern const double QR[80];
}  // testRiemannSolver

namespace testRiemannSolverLinear {
extern const double FL[80];
extern const double FR[80];
}  // namespace testRiemannSolverLinear

namespace testRiemannSolverNonlinear {
extern const double FL[80];
extern const double FR[80];
}  // namespace testRiemannNonlinear

namespace testSolutionUpdate {
extern const double luh[320];
}  // namespace testSolutionUpdate

namespace testSpaceTimePredictor {
extern const double luh[320];  // nVar * nDOFx * nDOFy * nDOFz
}  // namespace testSpaceTimePredictor

namespace testSpaceTimePredictorLinear {
extern const double lQhi[320];    // nVar * nDOFx * nDOFy * nDOFz
extern const double lQhbnd[480];  // nVar * nDOFy * nDOF_z * 6
}  // namespace testSpaceTimePredictorLinear

namespace testSpaceTimePredictorNonlinear {
extern const double lQi[1280];    // nVar * nDOFt * nDOFx * nDOFy * nDOFz
extern const double lFi[3840];    // nVar * nDOFx * nDOFy * nDOFz * nDOFt * dim
extern const double lQhi[320];    // nVar * nDOFx * nDOFy * nDOFz
extern const double lFhi[960];    // nVar * nDOFx * nDOFy * nDOFz * dim
extern const double lQhbnd[480];  // nVar * nDOFy * nDOF_z * 6
extern const double lFhbnd[480];  // nVar * nDOFy * nDOF_z * 6
}  // namespace testSpaceTimePredictorNonlinear

#endif  // Dim3

}  // namespace generic_euler
}  // namespace testdata
}  // namespace tests
}  // namespace exahype

#endif  // _EXAHYPE_TESTS_TESTDATA_GENERIC_EULER_TESTDATA_H_
