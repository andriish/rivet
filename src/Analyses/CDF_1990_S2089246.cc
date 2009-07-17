// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/CDF_1990_S2089246.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"

namespace Rivet {


  CDF_1990_S2089246::CDF_1990_S2089246()
    : Analysis("CDF_1990_S2089246")
  {
    setBeams(PROTON, ANTIPROTON);
    addProjection(ChargedFinalState(-3.5, 3.5), "FS");
    addProjection(Beam(), "Beam");
  }



  void CDF_1990_S2089246::init() {
    _hist_eta1800 = bookHistogram1D(3, 1, 1);
    _hist_eta630 = bookHistogram1D(4, 1, 1);
  }



  // Do the analysis
  void CDF_1990_S2089246::analyze(const Event& event) {
    const double sqrtS = applyProjection<Beam>(event, "Beam").sqrtS();
    const double weight = event.weight();

    const FinalState& fs = applyProjection<FinalState>(event, "FS");
    foreach (const Particle& p, fs.particles()) {
      const double eta = p.momentum().pseudorapidity();
      if (fuzzyEquals(sqrtS/GeV, 630)) {
        _hist_eta630->fill(eta, weight);
      } else if (fuzzyEquals(sqrtS/GeV, 1800)) {
        _hist_eta1800->fill(eta, weight);
      }
    }
  }

  
  
  // Finalize
  void CDF_1990_S2089246::finalize() {
    // Divide through by num events to get d<N>/d(eta) in bins
    scale(_hist_eta630, 1/sumOfWeights());
    scale(_hist_eta1800, 1/sumOfWeights());
  }
  
  
}
