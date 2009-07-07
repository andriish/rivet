// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/JetAlg.hh"

namespace Rivet {


  JetAlg::JetAlg(const FinalState& fs) {
    setName("JetAlg");
    VisibleFinalState vfs(fs);
    getLog() << Log::DEBUG << "Making visible final state from provided FS" << endl;
    addProjection(vfs, "FS");
  }
  
  
}
