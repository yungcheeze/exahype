#include "EulerFlow/tests/TestCase.h"


#include "tarch/compiler/CompilerSpecificSettings.h"
#include "tarch/tests/TestCaseFactory.h"
registerTest(exahype::tests::TestCase)


#ifdef UseTestSpecificCompilerSettings
#pragma optimize("",off)
#endif

 
exahype::tests::TestCase::TestCase():
  tarch::tests::TestCase( "exahype::tests::TestCase" ) {
}


exahype::tests::TestCase::~TestCase() {
}


void exahype::tests::TestCase::run() {
  // @todo If you have further tests, add them here
  testMethod( test1 );
  testMethod( test2 );
  testMethod( test3 );
}


void exahype::tests::TestCase::test1() {
  // @todo Add your test here
  validateEquals(1,1);



/*
  int trianglesA = 1;
  double xA[]    = {0.0, 0.0, 0.0};
  double yA[]    = {0.0, 1.0, 0.0};
  double zA[]    = {0.0, 0.0, 1.0};

  int trianglesB = 1;
  double xB[]    = {0.2, 0.2, 0.2};
  double yB[]    = {0.0, 1.0, 0.0};
  double zB[]    = {0.0, 0.0, 1.0};

  int numberOfContactPoints;

  double xC[]    = {-1.0, -1.0, -1.0};
  double yC[]    = {-1.0, -1.0, -1.0};
  double zC[]    = {-1.0, -1.0, -1.0};

  double xN[]    = {-1.0, -1.0, -1.0};
  double yN[]    = {-1.0, -1.0, -1.0};
  double zN[]    = {-1.0, -1.0, -1.0};

  ttd::ContactPoints::computeCollisionPoints(
    trianglesA, xA, yA, zA,
    trianglesB, xB, yB, zB,
    0.1,
    numberOfContactPoints,
    xC, yC, zC,
    xN, yN, zN
  );

  validateEquals(numberOfContactPoints,0);

  ttd::ContactPoints::computeCollisionPoints(
    trianglesA, xA, yA, zA,
    trianglesB, xB, yB, zB,
    0.25,
    numberOfContactPoints,
    xC, yC, zC,
    xN, yN, zN
  );

  validateEquals(numberOfContactPoints,1);

  validateNumericalEqualsWithParams3(xC[0],0.1,xC[0],yC[0],zC[0]);
  validateNumericalEqualsWithParams3(yC[0],0.0,xC[0],yC[0],zC[0]);
  validateNumericalEqualsWithParams3(zC[0],0.0,xC[0],yC[0],zC[0]);

  validateNumericalEqualsWithParams3(xN[0],0.1,xN[0],yN[0],zN[0]);
  validateNumericalEqualsWithParams3(yN[0],0.0,xN[0],yN[0],zN[0]);
  validateNumericalEqualsWithParams3(zN[0],0.0,xN[0],yN[0],zN[0]);
*/

}


void exahype::tests::TestCase::test2() {
  // @todo Add your test here
  validateEquals(2,2);
}


void exahype::tests::TestCase::test3() {
  // @todo Add your test here
  validateEquals(3,3);
}


#ifdef UseTestSpecificCompilerSettings
#pragma optimize("",on)
#endif
