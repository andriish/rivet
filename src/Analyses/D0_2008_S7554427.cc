// -*- C++ -*-
#include "Rivet/Analyses/D0_2008_S7554427.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  D0_2008_S7554427::D0_2008_S7554427()
  {
    // Run II Z pT
    setBeams(PROTON, ANTIPROTON);
    
    ZFinder zfinder(-MaxRapidity, MaxRapidity, 25.0*GeV, ELECTRON,
                    70.0*GeV, 110.0*GeV, 0.2);
    addProjection(zfinder, "ZFinder");
  } 



  // Book histograms
  void D0_2008_S7554427::init() {
    _h_ZpT         = bookHistogram1D(1, 1, 1, "Z pT",
                                     "$p_{\\perp}$(Z) [GeV]",
                                     "$1/\\sigma \\text{d}\\sigma/\\text{d}p_\\perp(Z)$");
    _h_forward_ZpT = bookHistogram1D(3, 1, 1, "Z pT (forward region only)",
                                     "$p_{\\perp}$(Z) [GeV]",
                                     "$1/\\sigma \\text{d}\\sigma/\\text{d}p_\\perp(Z)$");
  }



  // Do the analysis 
  void D0_2008_S7554427::analyze(const Event & e) {
    double weight = e.weight();

    const ZFinder& zfinder = applyProjection<ZFinder>(e, "ZFinder");
    if (zfinder.particles().size() == 1) {
      double yZ = fabs(zfinder.particles()[0].momentum().rapidity());
      double pTZ = zfinder.particles()[0].momentum().pT();
      _h_ZpT->fill(pTZ, weight);
      if (yZ > 2.0) {
        _h_forward_ZpT->fill(pTZ, weight);
      }
    }
    else {
      getLog() << Log::DEBUG << "no unique lepton pair found." << endl;
    }

  }



  // Finalize
  void D0_2008_S7554427::finalize() {
    normalize(_h_ZpT);
    normalize(_h_forward_ZpT);
  }

}
