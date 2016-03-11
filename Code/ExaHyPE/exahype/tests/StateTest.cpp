#include "exahype/tests/StateTest.h"

#include "tarch/compiler/CompilerSpecificSettings.h"
#include "tarch/tests/TestCaseFactory.h"
registerTest(exahype::tests::StateTest)
#ifdef UseTestSpecificCompilerSettings
#pragma optimize("", off)
#endif

    exahype::tests::StateTest::StateTest()
    : tarch::tests::TestCase("exahype::tests::StateTest") {
}

exahype::tests::StateTest::~StateTest() {}

void exahype::tests::StateTest::run() {
  // @todo If you have further tests, add them here
  testMethod(testState);
}

void exahype::tests::StateTest::testState() {
  // @todo Add your test here
  validateEquals(1, 1);
}

#ifdef UseTestSpecificCompilerSettings
#pragma optimize("", on)
#endif
