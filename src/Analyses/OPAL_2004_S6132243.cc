// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/OPAL_2004_S6132243.hh"


namespace Rivet {


  void OPAL_2004_S6132243::init() {
    // Book histograms
    // histChTot_       = bookHistogram1D("TotalChMult","Total charged multiplicity", 25, 1.0, 51.0);
    histChTot_       = bookHistogram1D(1, 1, 1, "Total charged multiplicity");
    histSphericity_  = bookHistogram1D("Sphericity", "Sphericity", 8, 0.0, 0.70);
    histAplanarity_  = bookHistogram1D("Aplanarity", "Aplanarity", 10, 0.0, 0.09);
    histPlanarity_   = bookHistogram1D("Planarity",  "Planarity", 16, 0.0, 0.70);
  }


  // Do the analysis
  void OPAL_2004_S6132243::analyze(const Event & event) {
    Log& log = getLog();

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
  }


  // Finalize
  void OPAL_2004_S6132243::finalize() { 
    // Normalize the histogram areas to 1
    normalize(histChTot_);
    normalize(histSphericity_);
    normalize(histPlanarity_); 
    normalize(histAplanarity_);
  }

}
