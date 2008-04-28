// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/ALEPH_1991_S2435284.hh"


namespace Rivet {

  void ALEPH_1991_S2435284::init() {
    // Book histogram
    _histChTot = bookHistogram1D(1, 1, 1, "Total charged multiplicity");
  }


  // Do the analysis
  void ALEPH_1991_S2435284::analyze(const Event& event) {
    Log& log = getLog();
    const Multiplicity& m = applyProjection<Multiplicity>(event, "Mult");
    log << Log::DEBUG << "Total charged multiplicity = " << m.totalMultiplicity() << endl;
    _histChTot->fill(m.totalMultiplicity(), event.weight());
  }


  // Finalize
  void ALEPH_1991_S2435284::finalize() { 
    // Normalize the histogram
    normalize(_histChTot);
  }

}
