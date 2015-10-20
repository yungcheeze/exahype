#include "ExplicitEulerForHeatEquation/tests/TestCase.h"


#include "tarch/compiler/CompilerSpecificSettings.h"
#include "tarch/tests/TestCaseFactory.h"
registerTest(myproject::tests::TestCase)


#ifdef UseTestSpecificCompilerSettings
#pragma optimize("",off)
#endif

 
myproject::tests::TestCase::TestCase():
  tarch::tests::TestCase( "myproject::tests::TestCase" ) {
}


myproject::tests::TestCase::~TestCase() {
}


void myproject::tests::TestCase::run() {
  // @todo If you have further tests, add them here
  testMethod( test1 );
  testMethod( test2 );
  testMethod( test3 );
}


void myproject::tests::TestCase::test1() {
  // @todo Add your test here
  validateEquals(1,1);
}


void myproject::tests::TestCase::test2() {
  // @todo Add your test here
  validateEquals(2,2);
}


void myproject::tests::TestCase::test3() {
  // @todo Add your test here
  validateEquals(3,3);
}


#ifdef UseTestSpecificCompilerSettings
#pragma optimize("",on)
#endif
