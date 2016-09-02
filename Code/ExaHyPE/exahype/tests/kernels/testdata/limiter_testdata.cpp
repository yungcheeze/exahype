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

#include "limiter_testdata.h"

namespace exahype {
namespace tests {
namespace testdata {
namespace limiter {
  
// nVar = 5
// basisSize = 4 (<=> N+1)
// basisSizeLim = 7 (<=> 2N+1)

#ifdef Dim2
namespace arraySize {
  const int sizeLuh = 80;
  const int sizeLim = 245;
}

namespace testFromLuhConversion {
  const double luh_in[80] = {};
  const double lob_out[80] = {}; 
  const double lim_out[245] = {}; 
}

namespace testToLuhConversion {
  const double lim_in[245] = {};
  const double luh_out[80] = {};
}

namespace testLocalMinMax {
  const double luh_in[80] = {};
  const double min_out[5] = {};
  const double max_out[5] = {};
}
#endif // Dim2

#ifdef Dim3
namespace arraySize {
  const int sizeLuh = 320;
  const int sizeLim = 1715;
}

namespace testFromLuhConversion {
  const double luh_in[320] = {};
  const double lob_out[320] = {};
  const double lim_out[1715] = {};
}

namespace testToLuhConversion {
  const double lim_in[1715] = {};
  const double luh_out[320] = {};
}

namespace testLocalMinMax {
  const double luh_in[320] = {};
  const double min_out[5] = {};
  const double max_out[5] = {};
}
#endif // Dim3


}  // namespace limiter
}  // namespace testdata
}  // namespace tests
}  // namespace exahype