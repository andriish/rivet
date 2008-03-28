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


  bool sortJetsByEt(const Jet& a, const Jet& b) {
    //for (Jet::const_iterator )
  }


  // Do the analysis
  void CDF_2004_S5839831::analyze(const Event& event) {
    Log log = getLog();

    const FastJets& jetproj = event.applyProjection(_jetproj);
    Jets jets = jetproj.getJets();
    
    /// @todo Sort jets by ET --- define ET!
    //sort(jets.begin(), jets.end(), sortJetsByEt);

    /// @todo Ensure there is only one well-defined primary vertex, i.e. no pileup. 
    /// Should be automatic as long as generator is run sensibly.

/*

    // Leading jet must be in central |eta| < 0.5 region.
    const double etaLead = jets.first().getPseudorapidity();
    if (fabs(etaLead) > 0.5) vetoEvent(event);
    
    // Get Et of the leading jet: used to bin histograms.
    const double ETlead = jets.first().getET();

    // Get the event weight.
    const double weight = event.weight();

    for (Jet j : jets) {

      const double ptMaxMin = ptMax + ptMin;
      if (inRange(etaLead/GeV, 40, 80) {
        _pt90Dbn1800Et40->fill(ptMaxMin/GeV, weight);
      else if (inRange(etaLead/GeV, 80, 120) {
        _pt90Dbn1800Et80->fill(ptMaxMin/GeV, weight);
      else if (inRange(etaLead/GeV, 120, 160) {
        _pt90Dbn1800Et120->fill(ptMaxMin/GeV, weight);
      else if (inRange(etaLead/GeV, 160, 200) {
        _pt90Dbn1800Et160->fill(ptMaxMin/GeV, weight);
      else if (inRange(etaLead/GeV, 200, 270) {
        _pt90Dbn1800Et200->fill(ptMaxMin/GeV, weight);
      }


    }

*/
  }


  void CDF_2004_S5839831::finalize() { 
    //const double avgNumTrans = _totalNumTrans / sumOfWeights();
    //normalize(_ptTrans2, avgNumTrans);
    //normalize(_ptTrans5, avgNumTrans);
    //normalize(_ptTrans30, avgNumTrans);
  }


}
