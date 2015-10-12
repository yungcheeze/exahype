#include "ExaHyPE/tests/TestCase.h"


#include "tarch/compiler/CompilerSpecificSettings.h"
#include "tarch/tests/TestCaseFactory.h"
registerTest(ExaHyPE::tests::TestCase)


#ifdef UseTestSpecificCompilerSettings
#pragma optimize("",off)
#endif

 
ExaHyPE::tests::TestCase::TestCase():
  tarch::tests::TestCase( "ExaHyPE::tests::TestCase" ) {
}


ExaHyPE::tests::TestCase::~TestCase() {
}


void ExaHyPE::tests::TestCase::run() {
  // @todo If you have further tests, add them here
  testMethod( test1 );
  testMethod( test2 );
  testMethod( test3 );
}


void ExaHyPE::tests::TestCase::test1() {
  // @todo Add your test here
  validateEquals(1,1);
}


void ExaHyPE::tests::TestCase::test2() {
  // @todo Add your test here
  validateEquals(2,2);
}


void ExaHyPE::tests::TestCase::test3() {
  // @todo Add your test here
  validateEquals(3,3);
}


#ifdef UseTestSpecificCompilerSettings
#pragma optimize("",on)
#endif
