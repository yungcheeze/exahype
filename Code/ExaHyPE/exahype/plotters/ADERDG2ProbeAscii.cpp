#include "exahype/plotters/ADERDG2ProbeAscii.h"

#include <fstream>
#include <cstring>

#include "peano/utils/Loop.h"

#include "kernels/DGBasisFunctions.h"


tarch::logging::Log exahype::plotters::ADERDG2ProbeAscii::_log( "exahype::plotters::ADERDG2ProbeAscii" );


exahype::plotters::ADERDG2ProbeAscii::ADERDG2ProbeAscii(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing):
  Device(postProcessing),
  _out(nullptr) {
}


exahype::plotters::ADERDG2ProbeAscii::~ADERDG2ProbeAscii() {
  if (_out!=nullptr) {
    _out->close();
    _out=nullptr;
  }
}


void exahype::plotters::ADERDG2ProbeAscii::startPlotting( double time ) {
  // In the very first time step, the
  if (time==std::numeric_limits<double>::max() ) {
    _time = 0.0;
  }
  else {
    _time = time;
  }
}


void exahype::plotters::ADERDG2ProbeAscii::finishPlotting() {
  if (_out!=nullptr && *_out) {
    (*_out) << std::endl;
  }
}


void exahype::plotters::ADERDG2ProbeAscii::init(const std::string& filename, int orderPlusOne, int unknowns, int writtenUnknowns, const std::string& select) {
  _order           = orderPlusOne-1;
  _solverUnknowns  = unknowns;
  _writtenUnknowns = writtenUnknowns;
  _select          = select;
  _filename        = filename;
  _time            = 0.0;

  _x(0) = Parser::getValueFromPropertyString( select, "x" );
  _x(1) = Parser::getValueFromPropertyString( select, "y" );
  #if DIM3
  _x(2) = Parser::getValueFromPropertyString( select, "z" );
  #endif

  logDebug( "init(...)", "probe at location " << _x );
}


void exahype::plotters::ADERDG2ProbeAscii::openOutputStream() {
  if (_out == nullptr) {
    if (tarch::la::oneEquals(_x,std::numeric_limits<double>::quiet_NaN())) {
      logError( "init(...)", "probe requires valid x, y (and z) coordinates in select statement. No plot written as plot location has been " << _x );
    }
    else {
      _out = new std::ofstream;

      std::ostringstream outputFilename;
      outputFilename << _filename
                   #ifdef Parallel
	           << "-rank-" << tarch::parallel::Node::getInstance().getRank()
                   #endif
		           << ".probe";
      _out->open( outputFilename.str() );
      
      // See issue #47 for discussion whether to quit program on failure
      if(*_out && _out->fail()) {
         logError("openOutputStream(...)", "Could not open file '" << outputFilename.str() << "': " << strerror(errno));
	 exit(-2);
      }

      if (*_out) {
        (*_out) << "# plot-time, real-time";
        for (int unknown=0; unknown < _writtenUnknowns; unknown++) {
          std::ostringstream identifier;
          identifier << "Q" << unknown;
          (*_out) << "," << identifier.str();
        }
        (*_out) << std::endl;
      }
    }
  }

  if (_out!=nullptr && *_out  ) {
    (*_out) << _time;
  }
}


std::string exahype::plotters::ADERDG2ProbeAscii::getIdentifier() {
  return "probe::ascii";
}


void exahype::plotters::ADERDG2ProbeAscii::plotPatch(
  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
  const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
  double timeStamp
) {
  if (
    tarch::la::allSmaller(offsetOfPatch,_x)
    &&
    tarch::la::allGreater(offsetOfPatch+sizeOfPatch,_x)
  ) {
    // lazy opening
    openOutputStream();

    // Map coordinate vector x onto reference element.
    tarch::la::Vector<DIMENSIONS,double> xRef = _x - offsetOfPatch;
    xRef(0) /=  sizeOfPatch(0);
    xRef(1) /=  sizeOfPatch(1);
    #ifdef Dim3
    xRef(2) /=  sizeOfPatch(2);
    #endif

    double* interpoland = new double[_solverUnknowns];
    double* value       = _writtenUnknowns==0 ? nullptr : new double[_writtenUnknowns];

    for (int unknown=0; unknown < _solverUnknowns; unknown++) {
      interpoland[unknown] = 0.0;

      // The code below evaluates the basis functions at the reference coordinates
      // and multiplies them with their respective coefficient.
      dfor(ii,_order+1) { // Gauss-Legendre node indices
        int iGauss = peano::utils::dLinearisedWithoutLookup(ii,_order + 1);
        interpoland[unknown] += kernels::basisFunctions[_order][ii(0)](xRef(0)) *
                 kernels::basisFunctions[_order][ii(1)](xRef(1)) *
                 #ifdef Dim3
                 kernels::basisFunctions[_order][ii(2)](xRef(2)) *
                 #endif
                 u[iGauss * _solverUnknowns + unknown];
        assertion3(interpoland[unknown] == interpoland[unknown], offsetOfPatch, sizeOfPatch, iGauss);
      }
    }

    _postProcessing->mapQuantities(
      offsetOfPatch,
      sizeOfPatch,
      _x,
      interpoland,
      value,
      timeStamp
    );

    if (_out!=nullptr && *_out) {
      (*_out) << ", " << timeStamp;
      for (int i=0; i<_writtenUnknowns; i++) {
        (*_out) << ", " << value[i];
      }
    }

    delete[] interpoland;
    delete[] value;
  }
}
