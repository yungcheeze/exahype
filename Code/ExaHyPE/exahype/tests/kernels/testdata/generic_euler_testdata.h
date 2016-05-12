#ifndef _EXAHYPE_TESTS_TESTDATA_GENERIC_EULER_TESTDATA_H_
#define _EXAHYPE_TESTS_TESTDATA_GENERIC_EULER_TESTDATA_H_

namespace exahype {
namespace tests {
namespace testdata {
namespace generic_euler {

#ifdef Dim3

namespace testPDEFluxes {
extern const double f[5];
extern const double g[5];
extern const double h[5];
}  // namespace testPDEFluxes

namespace testVolumeIntegral {
extern const double lduh[320];
}  // namespace testVolumeIntegral

namespace testSurfaceIntegral {
extern const double lduh[320];
}  // namesapce testSurfaceIntegral

namespace testSolutionUpdate {
extern const double luh[320];
}  // namespace testSolutionUpdate

#endif  // Dim3

}  // namespace generic_euler
}  // namespace testdata
}  // namespace tests
}  // namespace exahype

#endif  // _EXAHYPE_TESTS_TESTDATA_GENERIC_EULER_TESTDATA_H_
