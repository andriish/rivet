// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/PL273B181.hh"
#include "HepPDT/ParticleID.hh"

using namespace AIDA;
using namespace HepMC;


namespace Rivet {

  void PL273B181::init() {
    // Book histogram
    _histChTot = bookHistogram1D(1, 1, 1, "Total charged multiplicity");
  }


  // Do the analysis
  void PL273B181::analyze(const Event& event) {
    Log& log = getLog();
    const Multiplicity& m = event.applyProjection(_multproj);
    log << Log::INFO << "Total charged multiplicity = " << m.totalChargedMultiplicity() << endl;
    _histChTot->fill(m.totalChargedMultiplicity(), event.weight());
  }


  // Finalize
  void PL273B181::finalize() { 
    // Normalize the histogram
    normalize(_histChTot);
  }

}
