// -*- C++ -*-

#include "Rivet/Projections/Projection.hh"
#include "Rivet/Tools/Logging.hh"

using namespace Rivet;

Log& Projection::getLog() {
  string logname = "Rivet.Projection." + name();
  return Log::getLog(logname);
}
