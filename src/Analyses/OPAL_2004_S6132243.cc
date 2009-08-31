// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"

namespace Rivet {


  class OPAL_2004_S6132243 : public Analysis { 

    OPAL_2004_S6132243() : Analysis("OPAL_2004_S6132243") { }
    void init() { }
    void analyze(const Event & event) { }
    void finalize() { }

  };


  // This global object acts as a hook for the plugin system
  AnalysisBuilder<OPAL_2004_S6132243> plugin_OPAL_2004_S6132243;

}
