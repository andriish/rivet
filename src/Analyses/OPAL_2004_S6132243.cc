// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/OPAL_2004_S6132243.hh"

namespace Rivet {


  void OPAL_2004_S6132243::init() { }
  void OPAL_2004_S6132243::analyze(const Event & event) { }
  void OPAL_2004_S6132243::finalize() { }


  // This global object acts as a hook for the plugin system
  AnalysisBuilder<OPAL_2004_S6132243> plugin_OPAL_2004_S6132243;

}
