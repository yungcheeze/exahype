#include "exahype/plotters/ADERDG2AsciiSeismogram.h"


tarch::logging::Log exahype::plotters::ADERDG2AsciiSeismogram::_log( "exahype::plotters::ADERDG2AsciiSeismogram" );


exahype::plotters::ADERDG2AsciiSeismogram::ADERDG2AsciiSeismogram():
  _out(nullptr) {
}


exahype::plotters::ADERDG2AsciiSeismogram::~ADERDG2AsciiSeismogram() {
  if (_out!=nullptr) {
    _out->close();
    _out=nullptr;
  }
}


void exahype::plotters::ADERDG2AsciiSeismogram::startPlotting( double time ) {
  if (_out!=nullptr && *_out) {
    (*_out) << time;
  }
}


void exahype::plotters::ADERDG2AsciiSeismogram::finishPlotting() {
  if (_out!=nullptr && *_out) {
    (*_out) << std::endl;
  }
}


void exahype::plotters::ADERDG2AsciiSeismogram::init(const std::string& filename, int order, int unknowns, const std::string& select) {
  _order    = order;
  _unknowns = unknowns;
  _select   = select;

  _x(0) = getValueFromPropertyString( select, "x" );
  _x(1) = getValueFromPropertyString( select, "y" );
  #if DIM3
  _x(2) = getValueFromPropertyString( select, "z" );
  #endif

  logDebug( "init(...)", "probe at location " << _x );

  if ( !tarch::la::oneEquals(_x,std::numeric_limits<double>::quiet_NaN()) ) {
    _out = new std::ofstream;
    _out->open( filename );
  }

  if (_out!=nullptr && *_out) {
    (*_out) << "# time";
    for (int unknown=0; unknown < _unknowns; unknown++) {
      std::ostringstream identifier;
      identifier << "Q" << unknown;

      if ( _select.find(identifier.str())!=std::string::npos || _select.find("all")!=std::string::npos ) {
        (*_out) << "," << identifier.str();
      }
    }
    (*_out) << std::endl;
  }
}


std::string exahype::plotters::ADERDG2AsciiSeismogram::getIdentifier() {
  return "seismogram::ascii";
}


void exahype::plotters::ADERDG2AsciiSeismogram::plotPatch(
  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
  const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
  double timeStamp
) {

}
