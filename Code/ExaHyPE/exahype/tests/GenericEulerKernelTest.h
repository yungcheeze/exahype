// @todo ExaHyPE LIzenz
#ifndef _EXAHYPE_TESTS_GENERIC_EULER_KERNEL_TEST_H_
#define _EXAHYPE_TESTS_GENERIC_EULER_KERNEL_TEST_H_
 


#include "tarch/tests/TestCase.h"


namespace exahype {
    namespace tests {
      class GenericEulerKernelTest;
    } 
}
 

/**
 * This is just a default test case that demonstrated how to write unit tests 
 * in Peano. Feel free to rename, remove, or duplicate it. 
 */ 
class exahype::tests::GenericEulerKernelTest: public tarch::tests::TestCase {
  private:
    /**
     * These operation usually implement the real tests.
     */
    void testSpaceTimePredictor2d();

    const double eps = 1.0e-10; // for quick adaption of the test cases (say, switch to single precision)

    static void testFlux(const double* const Q, double* f, double* g);
  public: 
    GenericEulerKernelTest();
    virtual ~GenericEulerKernelTest();
     
    virtual void run();
};


#endif
