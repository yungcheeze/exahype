#include "exahype/tests/TestCase.h"

#include "tarch/compiler/CompilerSpecificSettings.h"
#include "tarch/tests/TestCaseFactory.h"
registerTest(exahype::tests::TestCase)
#ifdef UseTestSpecificCompilerSettings
#pragma optimize("", off)
#endif

    exahype::tests::TestCase::TestCase()
    : tarch::tests::TestCase("exahype::tests::TestCase") {
}

exahype::tests::TestCase::~TestCase() {}

void exahype::tests::TestCase::run() {
  // @todo If you have further tests, add them here
  testMethod(test1);
  testMethod(test2);
  testMethod(test3);
}

void exahype::tests::TestCase::test1() {
  // @todo Add your test here
  validateEquals(1, 1);
}

void exahype::tests::TestCase::test2() {
  // @todo Add your test here
  validateEquals(2, 2);
}

void exahype::tests::TestCase::test3() {
  // @todo Add your test here
  validateEquals(3, 3);
}

#ifdef UseTestSpecificCompilerSettings
#pragma optimize("", on)
#endif
