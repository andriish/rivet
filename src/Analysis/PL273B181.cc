// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analysis/PL273B181.hh"
#include "HepPDT/ParticleID.hh"

using namespace AIDA;
using namespace HepMC;

namespace Rivet {

  void PL273B181::init() {
    // Book histograms
    histChTot_       = bookHistogram1D("TotalChMult","Total charged multiplicity", 25, 1.0, 51.0);
    histSphericity_  = bookHistogram1D("Sphericity", "Event Shape: Sphericity", 8, 0.0, 0.70);
    histAplanarity_  = bookHistogram1D("Aplanarity", "Event Shape: APlanarity", 10, 0.0, 0.09);
    histPlanarity_   = bookHistogram1D("Planarity",  "Event Shape: Planarity", 16, 0.0, 0.70);
  }


  // Do the analysis
  void PL273B181::analyze(const Event & event) {
    Log& log = getLog();
    log << Log::DEBUG << "Starting analyzing" << endl;

    // Analyse and print some info
    const Multiplicity& m = event.applyProjection(mult);
    log << Log::INFO << "Total charged multiplicity    = " << m.totalChargedMultiplicity() << endl;

    //Analyse the event shape info
    const Sphericity& s = event.applyProjection(spher);
    log << Log::INFO << "Sphericity    = " << s.sphericity() << endl;
    log << Log::INFO << "Aplanarity    = " << s.aplanarity() << endl;
    log << Log::INFO << "Planarity     = " << s.planarity() << endl;

    // Fill histograms here, and scale them later
    const double weight = event.weight();
    histChTot_->fill(m.totalChargedMultiplicity(), weight);
    histSphericity_->fill(s.sphericity(), weight);
    histPlanarity_->fill(s.planarity(), weight);
    histAplanarity_->fill(s.aplanarity(), weight);

    // Finished...
    log << Log::DEBUG << "Finished analyzing" << endl;
  }


  // Finalize
  void PL273B181::finalize() { 
    // Normalize the histogram areas to 1
    normalize(histChTot_);
    normalize(histSphericity_);
    normalize(histPlanarity_); 
    normalize(histAplanarity_);
    //Log& log = getLog();
  }

}
