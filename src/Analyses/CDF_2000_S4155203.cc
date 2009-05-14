// -*- C++ -*-
// CDF Z pT analysis

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2000_S4155203.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  CDF_2000_S4155203::CDF_2000_S4155203() { 
    setBeams(PROTON, ANTIPROTON);
    ZFinder zfinder(FinalState(), ELECTRON, 66.0*GeV, 116.0*GeV, 0.2);
    addProjection(zfinder, "ZFinder");
  }


  // Book histograms
  void CDF_2000_S4155203::init() {
    /// @todo Cross-section units in label
    _hist_zpt =
      bookHistogram1D(1, 1, 1, "$p_\\perp$ of Z boson in $\\Pelectron \\Ppositron$ decays",
                      "$p_\\perp(\\PZ)$ / GeV", "$\\mathrm{d}{\\sigma}/\\mathrm{d}{p_\\perp(\\PZ)}$");
  }


  // Do the analysis
  void CDF_2000_S4155203::analyze(const Event& e) {
    const ZFinder& zfinder = applyProjection<ZFinder>(e, "ZFinder");
    if (zfinder.particles().size() != 1) {
      getLog() << Log::DEBUG << "No unique e+e- pair found" << endl;
      vetoEvent(e);
    }

    FourMomentum pZ = zfinder.particles()[0].momentum();    
    getLog() << Log::DEBUG << "Dilepton mass = " << pZ.mass()/GeV << " GeV"  << endl;
    getLog() << Log::DEBUG << "Dilepton pT   = " << pZ.pT()/GeV << " GeV" << endl;
    _hist_zpt->fill(pZ.pT()/GeV, e.weight());
  }


  void CDF_2000_S4155203::finalize() {
    // Normalize to the experimental cross-section
    /// @todo Get norm from generator cross-section
    normalize(_hist_zpt, 247.4);
  }


}
