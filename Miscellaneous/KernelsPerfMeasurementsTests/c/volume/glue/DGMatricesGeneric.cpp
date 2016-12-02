/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/

#include "matrices.h"

double*** KxiGeneric;

void freeDGMatricesGeneric() {
  const int MAX_ORDER = 9;

  // vectors
  for (int ii = 0; ii < MAX_ORDER + 1; ii++) {
    for (int jj = 0; jj < ii + 1; jj++) {
      delete[] KxiGeneric[ii][jj];
    }
    delete[] KxiGeneric[ii];  
  }
  delete[] KxiGeneric;
}

void initDGMatricesGeneric() {
  // \todo please implement!
  const int MAX_ORDER = 9;  // todo dangerous

  KxiGeneric = new double**[MAX_ORDER + 1];     // ***


  // vectors
  for (int ii = 0; ii < MAX_ORDER + 1; ii++) {
    KxiGeneric[ii] = new double*[ii + 1];
    for (int jj = 0; jj < ii + 1; jj++) {
      KxiGeneric[ii][jj] = new double[ii + 1];
    }
  }


  // Generated via Miscellaneous/aderdg/generateLookupTable.py!
  // N=0
  KxiGeneric[0][0][0] = 0.000000000000e+00;
  // N=1
  KxiGeneric[1][0][0] = -8.660254037844e-01;
  KxiGeneric[1][0][1] = -8.660254037844e-01;
  KxiGeneric[1][1][0] = 8.660254037844e-01;
  KxiGeneric[1][1][1] = 8.660254037844e-01;

  // N=2
  KxiGeneric[2][0][0] = -1.075828707280e+00;
  KxiGeneric[2][0][1] = -5.737753105492e-01;
  KxiGeneric[2][0][2] = 3.586095690933e-01;
  KxiGeneric[2][1][0] = 1.434438276373e+00;
  KxiGeneric[2][1][1] = 0.000000000000e+00;
  KxiGeneric[2][1][2] = -1.434438276373e+00;
  KxiGeneric[2][2][0] = -3.586095690933e-01;
  KxiGeneric[2][2][1] = 5.737753105492e-01;
  KxiGeneric[2][2][2] = 1.075828707280e+00;

  // N=4
  KxiGeneric[4][0][0] = -1.200518144782e+00;
  KxiGeneric[4][0][1] = -4.596060639300e-01;
  KxiGeneric[4][0][2] = 1.713312197563e-01;
  KxiGeneric[4][0][3] = -1.169847996191e-01;
  KxiGeneric[4][0][4] = 1.307284012760e-01;
  KxiGeneric[4][1][0] = 1.824799516392e+00;
  KxiGeneric[4][1][1] = -3.629695921020e-01;
  KxiGeneric[4][1][2] = -8.165764223060e-01;
  KxiGeneric[4][1][3] = 4.444344937741e-01;
  KxiGeneric[4][1][4] = -4.644712559813e-01;
  KxiGeneric[4][2][0] = -9.580242263155e-01;
  KxiGeneric[4][2][1] = 1.150025350187e+00;
  KxiGeneric[4][2][2] = -3.789561257387e-16;
  KxiGeneric[4][2][3] = -1.150025350187e+00;
  KxiGeneric[4][2][4] = 9.580242263155e-01;
  KxiGeneric[4][3][0] = 4.644712559813e-01;
  KxiGeneric[4][3][1] = -4.444344937741e-01;
  KxiGeneric[4][3][2] = 8.165764223060e-01;
  KxiGeneric[4][3][3] = 3.629695921020e-01;
  KxiGeneric[4][3][4] = -1.824799516392e+00;
  KxiGeneric[4][4][0] = -1.307284012760e-01;
  KxiGeneric[4][4][1] = 1.169847996191e-01;
  KxiGeneric[4][4][2] = -1.713312197563e-01;
  KxiGeneric[4][4][3] = 4.596060639300e-01;
  KxiGeneric[4][4][4] = 1.200518144782e+00;

  // N=5
  KxiGeneric[5][0][0] = -1.224169543734e+00;
  KxiGeneric[5][0][1] = -4.413288217193e-01;
  KxiGeneric[5][0][2] = 1.517970442050e-01;
  KxiGeneric[5][0][3] = -8.993719169484e-02;
  KxiGeneric[5][0][4] = 7.511859039061e-02;
  KxiGeneric[5][0][5] = -9.186600192800e-02;
  KxiGeneric[5][1][0] = 1.903292389724e+00;
  KxiGeneric[5][1][1] = -4.238415570051e-01;
  KxiGeneric[5][1][2] = -7.510722400319e-01;
  KxiGeneric[5][1][3] = 3.527291509737e-01;
  KxiGeneric[5][1][4] = -2.728043343245e-01;
  KxiGeneric[5][1][5] = 3.239594478789e-01;
  KxiGeneric[5][2][0] = -1.096959317636e+00;
  KxiGeneric[5][2][1] = 1.258536155414e+00;
  KxiGeneric[5][2][2] = -1.183945216212e-01;
  KxiGeneric[5][2][3] = -9.804616767269e-01;
  KxiGeneric[5][2][4] = 5.910515206232e-01;
  KxiGeneric[5][2][5] = -6.499299175974e-01;
  KxiGeneric[5][3][0] = 6.499299175974e-01;
  KxiGeneric[5][3][1] = -5.910515206232e-01;
  KxiGeneric[5][3][2] = 9.804616767269e-01;
  KxiGeneric[5][3][3] = 1.183945216212e-01;
  KxiGeneric[5][3][4] = -1.258536155414e+00;
  KxiGeneric[5][3][5] = 1.096959317636e+00;
  KxiGeneric[5][4][0] = -3.239594478789e-01;
  KxiGeneric[5][4][1] = 2.728043343245e-01;
  KxiGeneric[5][4][2] = -3.527291509737e-01;
  KxiGeneric[5][4][3] = 7.510722400319e-01;
  KxiGeneric[5][4][4] = 4.238415570051e-01;
  KxiGeneric[5][4][5] = -1.903292389724e+00;
  KxiGeneric[5][5][0] = 9.186600192800e-02;
  KxiGeneric[5][5][1] = -7.511859039061e-02;
  KxiGeneric[5][5][2] = 8.993719169484e-02;
  KxiGeneric[5][5][3] = -1.517970442050e-01;
  KxiGeneric[5][5][4] = 4.413288217193e-01;
  KxiGeneric[5][5][5] = 1.224169543734e+00;
  
  // N=6
  KxiGeneric[6][0][0] = -1.238935766264e+00;
  KxiGeneric[6][0][1] = -4.303826998789e-01;
  KxiGeneric[6][0][2] = 1.410455208324e-01;
  KxiGeneric[6][0][3] = -7.719765204350e-02;
  KxiGeneric[6][0][4] = 5.655161135900e-02;
  KxiGeneric[6][0][5] = -5.284240273195e-02;
  KxiGeneric[6][0][6] = 6.821403787966e-02;
  KxiGeneric[6][1][0] = 1.953026162578e+00;
  KxiGeneric[6][1][1] = -4.607770622774e-01;
  KxiGeneric[6][1][2] = -7.146663514283e-01;
  KxiGeneric[6][1][3] = 3.093553548663e-01;
  KxiGeneric[6][1][4] = -2.090887754710e-01;
  KxiGeneric[6][1][5] = 1.885998842134e-01;
  KxiGeneric[6][1][6] = -2.397926195872e-01;
  KxiGeneric[6][2][0] = -1.187709638233e+00;
  KxiGeneric[6][2][1] = 1.326175685572e+00;
  KxiGeneric[6][2][2] = -1.855211098817e-01;
  KxiGeneric[6][2][3] = -8.996222605947e-01;
  KxiGeneric[6][2][4] = 4.704134683009e-01;
  KxiGeneric[6][2][5] = -3.879970696836e-01;
  KxiGeneric[6][2][6] = 4.762072093626e-01;
  KxiGeneric[6][3][0] = 7.782478695739e-01;
  KxiGeneric[6][3][1] = -6.872555116179e-01;
  KxiGeneric[6][3][2] = 1.077018244667e+00;
  KxiGeneric[6][3][3] = -1.856111636271e-16;
  KxiGeneric[6][3][4] = -1.077018244667e+00;
  KxiGeneric[6][3][5] = 6.872555116179e-01;
  KxiGeneric[6][3][6] = -7.782478695739e-01;
  KxiGeneric[6][4][0] = -4.762072093626e-01;
  KxiGeneric[6][4][1] = 3.879970696836e-01;
  KxiGeneric[6][4][2] = -4.704134683009e-01;
  KxiGeneric[6][4][3] = 8.996222605947e-01;
  KxiGeneric[6][4][4] = 1.855211098817e-01;
  KxiGeneric[6][4][5] = -1.326175685572e+00;
  KxiGeneric[6][4][6] = 1.187709638233e+00;
  KxiGeneric[6][5][0] = 2.397926195872e-01;
  KxiGeneric[6][5][1] = -1.885998842134e-01;
  KxiGeneric[6][5][2] = 2.090887754710e-01;
  KxiGeneric[6][5][3] = -3.093553548663e-01;
  KxiGeneric[6][5][4] = 7.146663514283e-01;
  KxiGeneric[6][5][5] = 4.607770622774e-01;
  KxiGeneric[6][5][6] = -1.953026162578e+00;
  KxiGeneric[6][6][0] = -6.821403787966e-02;
  KxiGeneric[6][6][1] = 5.284240273195e-02;
  KxiGeneric[6][6][2] = -5.655161135900e-02;
  KxiGeneric[6][6][3] = 7.719765204350e-02;
  KxiGeneric[6][6][4] = -1.410455208324e-01;
  KxiGeneric[6][6][5] = 4.303826998789e-01;
  KxiGeneric[6][6][6] = 1.238935766264e+00;

  // N=7
  KxiGeneric[7][0][0] = -1.248773141877e+00;
  KxiGeneric[7][0][1] = -4.232798984656e-01;
  KxiGeneric[7][0][2] = 1.344197824361e-01;
  KxiGeneric[7][0][3] = -7.000333367538e-02;
  KxiGeneric[7][0][4] = 4.754856157327e-02;
  KxiGeneric[7][0][5] = -3.933175770947e-02;
  KxiGeneric[7][0][6] = 3.941958371330e-02;
  KxiGeneric[7][0][7] = -5.270728187195e-02;
  KxiGeneric[7][1][0] = 1.986471526627e+00;
  KxiGeneric[7][1][1] = -4.849509784821e-01;
  KxiGeneric[7][1][2] = -6.920713528237e-01;
  KxiGeneric[7][1][3] = 2.847468046774e-01;
  KxiGeneric[7][1][4] = -1.781610100961e-01;
  KxiGeneric[7][1][5] = 1.419182262877e-01;
  KxiGeneric[7][1][6] = -1.395697200511e-01;
  KxiGeneric[7][1][7] = 1.849978723814e-01;
  KxiGeneric[7][2][0] = -1.249883197126e+00;
  KxiGeneric[7][2][1] = 1.371205081828e+00;
  KxiGeneric[7][2][2] = -2.277693241176e-01;
  KxiGeneric[7][2][3] = -8.533399324597e-01;
  KxiGeneric[7][2][4] = 4.117619919744e-01;
  KxiGeneric[7][2][5] = -2.984655560328e-01;
  KxiGeneric[7][2][6] = 2.811834246508e-01;
  KxiGeneric[7][2][7] = -3.657207457383e-01;
  KxiGeneric[7][3][0] = 8.690248642539e-01;
  KxiGeneric[7][3][1] = -7.532127092141e-01;
  KxiGeneric[7][3][2] = 1.139276942272e+00;
  KxiGeneric[7][3][3] = -6.884529508726e-02;
  KxiGeneric[7][3][4] = -9.885912999966e-01;
  KxiGeneric[7][3][5] = 5.497351352211e-01;
  KxiGeneric[7][3][6] = -4.712717926469e-01;
  KxiGeneric[7][3][7] = 5.902702071060e-01;
  KxiGeneric[7][4][0] = -5.902702071060e-01;
  KxiGeneric[7][4][1] = 4.712717926469e-01;
  KxiGeneric[7][4][2] = -5.497351352211e-01;
  KxiGeneric[7][4][3] = 9.885912999966e-01;
  KxiGeneric[7][4][4] = 6.884529508726e-02;
  KxiGeneric[7][4][5] = -1.139276942272e+00;
  KxiGeneric[7][4][6] = 7.532127092141e-01;
  KxiGeneric[7][4][7] = -8.690248642539e-01;
  KxiGeneric[7][5][0] = 3.657207457383e-01;
  KxiGeneric[7][5][1] = -2.811834246508e-01;
  KxiGeneric[7][5][2] = 2.984655560328e-01;
  KxiGeneric[7][5][3] = -4.117619919744e-01;
  KxiGeneric[7][5][4] = 8.533399324597e-01;
  KxiGeneric[7][5][5] = 2.277693241176e-01;
  KxiGeneric[7][5][6] = -1.371205081828e+00;
  KxiGeneric[7][5][7] = 1.249883197126e+00;
  KxiGeneric[7][6][0] = -1.849978723814e-01;
  KxiGeneric[7][6][1] = 1.395697200511e-01;
  KxiGeneric[7][6][2] = -1.419182262877e-01;
  KxiGeneric[7][6][3] = 1.781610100961e-01;
  KxiGeneric[7][6][4] = -2.847468046774e-01;
  KxiGeneric[7][6][5] = 6.920713528237e-01;
  KxiGeneric[7][6][6] = 4.849509784821e-01;
  KxiGeneric[7][6][7] = -1.986471526627e+00;
  KxiGeneric[7][7][0] = 5.270728187195e-02;
  KxiGeneric[7][7][1] = -3.941958371330e-02;
  KxiGeneric[7][7][2] = 3.933175770947e-02;
  KxiGeneric[7][7][3] = -4.754856157327e-02;
  KxiGeneric[7][7][4] = 7.000333367538e-02;
  KxiGeneric[7][7][5] = -1.344197824361e-01;
  KxiGeneric[7][7][6] = 4.232798984656e-01;
  KxiGeneric[7][7][7] = 1.248773141877e+00;

  // N=8
  KxiGeneric[8][0][0] = -1.255656088098e+00;
  KxiGeneric[8][0][1] = -4.183978248600e-01;
  KxiGeneric[8][0][2] = 1.300181658087e-01;
  KxiGeneric[8][0][3] = -6.548023049355e-02;
  KxiGeneric[8][0][4] = 4.236031185339e-02;
  KxiGeneric[8][0][5] = -3.262358480478e-02;
  KxiGeneric[8][0][6] = 2.916728809238e-02;
  KxiGeneric[8][0][7] = -3.064117428684e-02;
  KxiGeneric[8][0][8] = 4.197362432633e-02;
  KxiGeneric[8][1][0] = 2.010021351292e+00;
  KxiGeneric[8][1][1] = -5.016657854850e-01;
  KxiGeneric[8][1][2] = -6.769832419076e-01;
  KxiGeneric[8][1][3] = 2.692138635353e-01;
  KxiGeneric[8][1][4] = -1.602985146075e-01;
  KxiGeneric[8][1][5] = 1.187447077945e-01;
  KxiGeneric[8][1][6] = -1.039993129635e-01;
  KxiGeneric[8][1][7] = 1.080391382041e-01;
  KxiGeneric[8][1][8] = -1.472029988823e-01;
  KxiGeneric[8][2][0] = -1.294202075032e+00;
  KxiGeneric[8][2][1] = 1.402699759000e+00;
  KxiGeneric[8][2][2] = -2.562639281818e-01;
  KxiGeneric[8][2][3] = -8.239035397637e-01;
  KxiGeneric[8][2][4] = 3.777470995127e-01;
  KxiGeneric[8][2][5] = -2.540518734313e-01;
  KxiGeneric[8][2][6] = 2.124411755333e-01;
  KxiGeneric[8][2][7] = -2.154851142534e-01;
  KxiGeneric[8][2][8] = 2.903314666640e-01;
  KxiGeneric[8][3][0] = 9.350501663231e-01;
  KxiGeneric[8][3][1] = -8.002227249223e-01;
  KxiGeneric[8][3][2] = 1.181960666698e+00;
  KxiGeneric[8][3][3] = -1.131793140867e-01;
  KxiGeneric[8][3][4] = -9.369707771701e-01;
  KxiGeneric[8][3][5] = 4.816403690687e-01;
  KxiGeneric[8][3][6] = -3.644593173890e-01;
  KxiGeneric[8][3][7] = 3.529618140521e-01;
  KxiGeneric[8][3][8] = -4.658610418418e-01;
  KxiGeneric[8][4][0] = -6.759723042189e-01;
  KxiGeneric[8][4][1] = 5.324612399834e-01;
  KxiGeneric[8][4][2] = -6.055818291441e-01;
  KxiGeneric[8][4][3] = 1.047058839436e+00;
  KxiGeneric[8][4][4] = 4.399672026717e-16;
  KxiGeneric[8][4][5] = -1.047058839436e+00;
  KxiGeneric[8][4][6] = 6.055818291441e-01;
  KxiGeneric[8][4][7] = -5.324612399834e-01;
  KxiGeneric[8][4][8] = 6.759723042189e-01;
  KxiGeneric[8][5][0] = 4.658610418418e-01;
  KxiGeneric[8][5][1] = -3.529618140521e-01;
  KxiGeneric[8][5][2] = 3.644593173890e-01;
  KxiGeneric[8][5][3] = -4.816403690687e-01;
  KxiGeneric[8][5][4] = 9.369707771701e-01;
  KxiGeneric[8][5][5] = 1.131793140867e-01;
  KxiGeneric[8][5][6] = -1.181960666698e+00;
  KxiGeneric[8][5][7] = 8.002227249223e-01;
  KxiGeneric[8][5][8] = -9.350501663231e-01;
  KxiGeneric[8][6][0] = -2.903314666640e-01;
  KxiGeneric[8][6][1] = 2.154851142534e-01;
  KxiGeneric[8][6][2] = -2.124411755333e-01;
  KxiGeneric[8][6][3] = 2.540518734313e-01;
  KxiGeneric[8][6][4] = -3.777470995127e-01;
  KxiGeneric[8][6][5] = 8.239035397637e-01;
  KxiGeneric[8][6][6] = 2.562639281818e-01;
  KxiGeneric[8][6][7] = -1.402699759000e+00;
  KxiGeneric[8][6][8] = 1.294202075032e+00;
  KxiGeneric[8][7][0] = 1.472029988823e-01;
  KxiGeneric[8][7][1] = -1.080391382041e-01;
  KxiGeneric[8][7][2] = 1.039993129635e-01;
  KxiGeneric[8][7][3] = -1.187447077945e-01;
  KxiGeneric[8][7][4] = 1.602985146075e-01;
  KxiGeneric[8][7][5] = -2.692138635353e-01;
  KxiGeneric[8][7][6] = 6.769832419076e-01;
  KxiGeneric[8][7][7] = 5.016657854850e-01;
  KxiGeneric[8][7][8] = -2.010021351292e+00;
  KxiGeneric[8][8][0] = -4.197362432633e-02;
  KxiGeneric[8][8][1] = 3.064117428684e-02;
  KxiGeneric[8][8][2] = -2.916728809238e-02;
  KxiGeneric[8][8][3] = 3.262358480478e-02;
  KxiGeneric[8][8][4] = -4.236031185339e-02;
  KxiGeneric[8][8][5] = 6.548023049355e-02;
  KxiGeneric[8][8][6] = -1.300181658087e-01;
  KxiGeneric[8][8][7] = 4.183978248600e-01;
  KxiGeneric[8][8][8] = 1.255656088098e+00;

  //N=9
  KxiGeneric[9][0][0] = -1.260660205839e+00;
  KxiGeneric[9][0][1] = -4.148927713020e-01;
  KxiGeneric[9][0][2] = 1.269322840112e-01;
  KxiGeneric[9][0][3] = -6.242556341322e-02;
  KxiGeneric[9][0][4] = 3.904733259262e-02;
  KxiGeneric[9][0][5] = -2.869242542062e-02;
  KxiGeneric[9][0][6] = 2.397617139031e-02;
  KxiGeneric[9][0][7] = -2.260981542902e-02;
  KxiGeneric[9][0][8] = 2.455626988016e-02;
  KxiGeneric[9][0][9] = -3.422882091683e-02;
  KxiGeneric[9][1][0] = 2.027220013968e+00;
  KxiGeneric[9][1][1] = -5.137174273279e-01;
  KxiGeneric[9][1][2] = -6.663645647310e-01;
  KxiGeneric[9][1][3] = 2.586902796970e-01;
  KxiGeneric[9][1][4] = -1.488666530544e-01;
  KxiGeneric[9][1][5] = 1.051510984479e-01;
  KxiGeneric[9][1][6] = -8.600065857106e-02;
  KxiGeneric[9][1][7] = 8.010053769945e-02;
  KxiGeneric[9][1][8] = -8.638173508757e-02;
  KxiGeneric[9][1][9] = 1.199851267913e-01;
  KxiGeneric[9][2][0] = -1.326846946353e+00;
  KxiGeneric[9][2][1] = 1.425593189709e+00;
  KxiGeneric[9][2][2] = -2.764647873771e-01;
  KxiGeneric[9][2][3] = -8.038396407681e-01;
  KxiGeneric[9][2][4] = 3.558856451654e-01;
  KxiGeneric[9][2][5] = -2.279530855410e-01;
  KxiGeneric[9][2][6] = 1.777094388193e-01;
  KxiGeneric[9][2][7] = -1.612329092336e-01;
  KxiGeneric[9][2][8] = 1.713638255696e-01;
  KxiGeneric[9][2][9] = -2.363446367747e-01;
  KxiGeneric[9][3][0] = 9.843513264005e-01;
  KxiGeneric[9][3][1] = -8.348390467020e-01;
  KxiGeneric[9][3][2] = 1.212574259881e+00;
  KxiGeneric[9][3][3] = -1.436881209529e-01;
  KxiGeneric[9][3][4] = -9.035733595898e-01;
  KxiGeneric[9][3][5] = 4.415232856090e-01;
  KxiGeneric[9][3][6] = -3.106478783087e-01;
  KxiGeneric[9][3][7] = 2.680707448618e-01;
  KxiGeneric[9][3][8] = -2.775392562152e-01;
  KxiGeneric[9][3][9] = 3.780658887103e-01;
  KxiGeneric[9][4][0] = -7.413087738687e-01;
  KxiGeneric[9][4][1] = 5.784158377192e-01;
  KxiGeneric[9][4][2] = -6.463526073133e-01;
  KxiGeneric[9][4][3] = 1.087886347304e+00;
  KxiGeneric[9][4][4] = -4.499318311735e-02;
  KxiGeneric[9][4][5] = -9.925290910988e-01;
  KxiGeneric[9][4][6] = 5.315862285369e-01;
  KxiGeneric[9][4][7] = -4.140039734282e-01;
  KxiGeneric[9][4][8] = 4.085606779488e-01;
  KxiGeneric[9][4][9] = -5.447221435019e-01;
  KxiGeneric[9][5][0] = 5.447221435019e-01;
  KxiGeneric[9][5][1] = -4.085606779488e-01;
  KxiGeneric[9][5][2] = 4.140039734282e-01;
  KxiGeneric[9][5][3] = -5.315862285369e-01;
  KxiGeneric[9][5][4] = 9.925290910988e-01;
  KxiGeneric[9][5][5] = 4.499318311735e-02;
  KxiGeneric[9][5][6] = -1.087886347304e+00;
  KxiGeneric[9][5][7] = 6.463526073133e-01;
  KxiGeneric[9][5][8] = -5.784158377192e-01;
  KxiGeneric[9][5][9] = 7.413087738687e-01;
  KxiGeneric[9][6][0] = -3.780658887103e-01;
  KxiGeneric[9][6][1] = 2.775392562152e-01;
  KxiGeneric[9][6][2] = -2.680707448618e-01;
  KxiGeneric[9][6][3] = 3.106478783087e-01;
  KxiGeneric[9][6][4] = -4.415232856090e-01;
  KxiGeneric[9][6][5] = 9.035733595898e-01;
  KxiGeneric[9][6][6] = 1.436881209529e-01;
  KxiGeneric[9][6][7] = -1.212574259881e+00;
  KxiGeneric[9][6][8] = 8.348390467020e-01;
  KxiGeneric[9][6][9] = -9.843513264005e-01;
  KxiGeneric[9][7][0] = 2.363446367747e-01;
  KxiGeneric[9][7][1] = -1.713638255696e-01;
  KxiGeneric[9][7][2] = 1.612329092336e-01;
  KxiGeneric[9][7][3] = -1.777094388193e-01;
  KxiGeneric[9][7][4] = 2.279530855410e-01;
  KxiGeneric[9][7][5] = -3.558856451654e-01;
  KxiGeneric[9][7][6] = 8.038396407681e-01;
  KxiGeneric[9][7][7] = 2.764647873771e-01;
  KxiGeneric[9][7][8] = -1.425593189709e+00;
  KxiGeneric[9][7][9] = 1.326846946353e+00;
  KxiGeneric[9][8][0] = -1.199851267913e-01;
  KxiGeneric[9][8][1] = 8.638173508757e-02;
  KxiGeneric[9][8][2] = -8.010053769945e-02;
  KxiGeneric[9][8][3] = 8.600065857106e-02;
  KxiGeneric[9][8][4] = -1.051510984479e-01;
  KxiGeneric[9][8][5] = 1.488666530544e-01;
  KxiGeneric[9][8][6] = -2.586902796970e-01;
  KxiGeneric[9][8][7] = 6.663645647310e-01;
  KxiGeneric[9][8][8] = 5.137174273279e-01;
  KxiGeneric[9][8][9] = -2.027220013968e+00;
  KxiGeneric[9][9][0] = 3.422882091683e-02;
  KxiGeneric[9][9][1] = -2.455626988016e-02;
  KxiGeneric[9][9][2] = 2.260981542902e-02;
  KxiGeneric[9][9][3] = -2.397617139031e-02;
  KxiGeneric[9][9][4] = 2.869242542062e-02;
  KxiGeneric[9][9][5] = -3.904733259262e-02;
  KxiGeneric[9][9][6] = 6.242556341322e-02;
  KxiGeneric[9][9][7] = -1.269322840112e-01;
  KxiGeneric[9][9][8] = 4.148927713020e-01;
  KxiGeneric[9][9][9] = 1.260660205839e+00;

}
