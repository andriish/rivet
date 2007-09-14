// -*- C++ -*-
#include "Rivet/Projection.hh"
#include "Rivet/Tools/Logging.hh"

using namespace Rivet;

Log& Projection::getLog() {
  string logname = "Rivet.Projection." + getName();
  return Log::getLog(logname);
}
