#include "EulerFlow3d/dg/DGHelpers.h"

const double exahype::dg::normal[3][6] = { { -1, 1,  0, 0,  0, 0},
                                           {  0, 0, -1, 1,  0, 0},
                                           {  0, 0,  0, 0, -1, 1} };
