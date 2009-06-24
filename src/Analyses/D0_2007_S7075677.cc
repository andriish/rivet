// -*- C++ -*-
#include "Rivet/Analyses/D0_2007_S7075677.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  D0_2007_S7075677::D0_2007_S7075677()
  {
    // Run II Z rapidity
    setBeams(PROTON, ANTIPROTON);
    
    std::vector<std::pair<double, double> > etaRanges;
    // remove eta cuts for the moment, because it seems, like they have been
    // corrected for.
    // todo: ask gavin hesketh about it, he first implemented this
    // analysis without eta cuts.
    //    etaRanges.push_back(make_pair(-3.2, -1.5));
    //    etaRanges.push_back(make_pair(-0.9, 0.9));
    //    etaRanges.push_back(make_pair(1.5, 3.2));
    ZFinder zfinder(etaRanges, 15.0*GeV, ELECTRON, 71.0*GeV, 111.0*GeV, 0.2);
    addProjection(zfinder, "ZFinder");
  } 



  // Book histograms
  void D0_2007_S7075677::init() {
    _h_yZ = bookHistogram1D(1, 1, 1, "Inclusive Z boson rapidity",
                            "$|y|$(Z)", "$1/\\sigma \\; \\text{d}\\sigma/\\text{d}|y|(Z)$");
  }



  // Do the analysis 
  void D0_2007_S7075677::analyze(const Event & e) {
    double weight = e.weight();

    const ZFinder& zfinder = applyProjection<ZFinder>(e, "ZFinder");
    if (zfinder.particles().size() == 1) {
      const ParticleVector& el(zfinder.constituentsFinalState().particles());
      if (el[0].momentum().pT() > 25.0*GeV || el[1].momentum().pT() > 25.0*GeV) {
        double yZ = fabs(zfinder.particles()[0].momentum().rapidity());
        _h_yZ->fill(yZ, weight);
      }
    }
    else {
      getLog() << Log::DEBUG << "no unique lepton pair found." << endl;
    }
  }



  // Finalize
  void D0_2007_S7075677::finalize() {
    // Data seems to have been normalized for the avg of the two sides 
    // (+ve & -ve rapidity) rather than the sum, hence the 0.5:
    normalize(_h_yZ, 0.5);
  }


}
