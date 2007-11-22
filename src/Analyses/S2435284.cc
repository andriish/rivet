// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/S2435284.hh"
#include "HepPDT/ParticleID.hh"

using namespace AIDA;
using namespace HepMC;


namespace Rivet {

  void S2435284::init() {
    // Book histogram
    _histChTot = bookHistogram1D(1, 1, 1, "Total charged multiplicity");
  }


  // Do the analysis
  void S2435284::analyze(const Event& event) {
    Log& log = getLog();
    const Multiplicity& m = event.applyProjection(_multproj);
    log << Log::DEBUG << "Total charged multiplicity = " << m.totalChargedMultiplicity() << endl;
    _histChTot->fill(m.totalChargedMultiplicity(), event.weight());
  }


  // Finalize
  void S2435284::finalize() { 
    // Normalize the histogram
    normalize(_histChTot);
  }

}
