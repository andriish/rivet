// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Projection class.
//

#include "Rivet/Projections/Projection.hh"
#include "Rivet/Tools/Logging.hh"

using namespace Rivet;

Projection::~Projection() {}

RivetInfo Projection::getInfo() const {
  return info;
}

Log& Projection::getLog() {
  string logname = "Rivet.Projection." + name();
  return Log::getLog(logname);
}
