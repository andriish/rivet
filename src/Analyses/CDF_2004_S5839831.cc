// -*- C++ -*-
// Underlying event analysis at CDF.

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2004_S5839831.hh"
#include "Rivet/RivetAIDA.hh"


namespace Rivet {


  // Book histograms
  void CDF_2004_S5839831::init() {
    _pt90Max1800 = bookProfile1D(2, 1, 1, "pTmax vs ET at sqrt{s} = 1800 GeV");
    _pt90Min1800 = bookProfile1D(2, 1, 2, "pTmin vs ET at sqrt{s} = 1800 GeV");
    _pt90Diff1800 = bookProfile1D(2, 1, 3, "pTdiff vs ET at sqrt{s} = 1800 GeV");
    _pt90Dbn1800Et40 = bookHistogram1D(3, 1, 1, "pT distribution in MAX+MIN transverse cones for 40 < ET < 80 GeV at sqrt{s} = 1800 GeV");
    _pt90Dbn1800Et80 = bookHistogram1D(3, 1, 2, "pT distribution in MAX+MIN transverse cones for 80 < ET < 120 GeV at sqrt{s} = 1800 GeV");
    _pt90Dbn1800Et120 = bookHistogram1D(3, 1, 3, "pT distribution in MAX+MIN transverse cones for 120 < ET < 160 GeV at sqrt{s} = 1800 GeV");
    _pt90Dbn1800Et160 = bookHistogram1D(3, 1, 4, "pT distribution in MAX+MIN transverse cones for 160 < ET < 200 GeV at sqrt{s} = 1800 GeV");
    _pt90Dbn1800Et200 = bookHistogram1D(3, 1, 5, "pT distribution in MAX+MIN transverse cones for 200 < ET < 270 GeV at sqrt{s} = 1800 GeV");
    _num90Max1800 = bookProfile1D(4, 1, 1, "Nmax vs ET at sqrt{s} = 1800 GeV");
    _num90Min1800 = bookProfile1D(4, 1, 2, "Nmin vs ET at sqrt{s} = 1800 GeV");    
    _numTracksDbn1800 = bookHistogram1D(5, 1, 1, "Track multiplicity distribution at sqrt{s} = 1800 GeV");
    _ptDbn1800 = bookHistogram1D(6, 1, 1, "pT distribution at sqrt{s} = 1800 GeV");
    _pTSum1800_2Jet = bookProfile1D(7, 1, 1, "pTsum vs ET (for removal of 2 jets) at sqrt{s} = 1800 GeV");
    _pTSum1800_3Jet = bookProfile1D(7, 1, 2, "pTsum vs ET (for removal of 3 jets) at sqrt{s} = 1800 GeV");            
    _pt90Max630 = bookProfile1D(8, 1, 1, "pTmax vs ET at sqrt{s} = 630 GeV");
    _pt90Min630 = bookProfile1D(8, 1, 2, "pTmin vs ET at sqrt{s} = 630 GeV");
    _pt90Diff630 = bookProfile1D(8, 1, 3, "pTdiff vs ET at sqrt{s} = 630 GeV");
    _pTSum630_2Jet = bookProfile1D(9, 1, 1, "pTsum vs ET (for removal of 2 jets) at sqrt{s} = 630 GeV");
    _pTSum630_3Jet = bookProfile1D(9, 1, 2, "pTsum vs ET (for removal of 3 jets) at sqrt{s} = 630 GeV");
  }


  // Do the analysis
  void CDF_2004_S5839831::analyze(const Event& event) {
    //Log log = getLog();

  }


  void CDF_2004_S5839831::finalize() { 
    //const double avgNumTrans = _totalNumTrans / sumOfWeights();
    //normalize(_ptTrans2, avgNumTrans);
    //normalize(_ptTrans5, avgNumTrans);
    //normalize(_ptTrans30, avgNumTrans);
  }


}
