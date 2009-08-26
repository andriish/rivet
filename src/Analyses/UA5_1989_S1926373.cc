// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/UA5_1989_S1926373.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {



  UA5_1989_S1926373::UA5_1989_S1926373()
    : Analysis("UA5_1989_S1926373")
  { 
    setBeams(PROTON, ANTIPROTON);
    addProjection(Beam(), "Beams");
    addProjection(ChargedFinalState(), "CFSAll");
    addProjection(ChargedFinalState(-0.5, 0.5), "CFS05");
    addProjection(ChargedFinalState(-1.5, 1.5), "CFS15");
    addProjection(ChargedFinalState(-3.0, 3.0), "CFS30");
    addProjection(ChargedFinalState(-5.0, 5.0), "CFS50");
  }



  // Book histograms
  void UA5_1989_S1926373::init() {
      // NB. _hist_nch{200,900} and _hist_nch{200,900}eta50 use the same data but different binning
     _hist_nch200       = bookHistogram1D(1,1,1); 
     _hist_nch900       = bookHistogram1D(2,1,1);
     _hist_nch200eta05  = bookHistogram1D(3,1,1);
     _hist_nch200eta15  = bookHistogram1D(4,1,1);
     _hist_nch200eta30  = bookHistogram1D(5,1,1);
     _hist_nch200eta50  = bookHistogram1D(6,1,1);
     _hist_nch900eta05  = bookHistogram1D(7,1,1);
     _hist_nch900eta15  = bookHistogram1D(8,1,1);
     _hist_nch900eta30  = bookHistogram1D(9,1,1);
     _hist_nch900eta50  = bookHistogram1D(10,1,1);
     _hist_mean_nch_200 = bookHistogram1D(11,1,1); 
     _hist_mean_nch_900 = bookHistogram1D(12,1,1);
  } 



  // Do the analysis
  void UA5_1989_S1926373::analyze(const Event& event) {
    const double sqrtS = applyProjection<Beam>(event, "Beams").sqrtS();
    const double weight = event.weight();
    
    // Minimum Bias trigger requirements from the hodoscopes
    int n_trig_1(0), n_trig_2(0);
    /// @todo Use CFS in +,- eta ranges as below, to cache this loop between UA5 analyses
    const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(event, "CFSAll");
    foreach (const Particle& p, cfs.particles()) {
      const double eta = p.momentum().pseudorapidity();
      if (inRange(eta, -5.6, -2.0)) n_trig_1 += 1;
      else if (inRange(eta, 2.0, 5.6)) n_trig_2 += 1;
    }
    
    // Require at least one coincidence hit in trigger hodoscopes
    getLog() << Log::DEBUG << "Trigger -: " << n_trig_1 << ", Trigger +: " << n_trig_2 << endl;
    if (!n_trig_1 || !n_trig_2) vetoEvent;

    // Count  final state particles in several eta regions
    const int numP05 = applyProjection<ChargedFinalState>(event, "CFS05").size();
    const int numP15 = applyProjection<ChargedFinalState>(event, "CFS15").size();
    const int numP30 = applyProjection<ChargedFinalState>(event, "CFS35").size();
    const int numP50 = applyProjection<ChargedFinalState>(event, "CFS50").size();
        
    // Fill histograms
    if (fuzzyEquals(sqrtS, 200.0, 1E-4)) {
      _hist_nch200->fill(numP50, weight);
      _hist_nch200eta05->fill(numP05, weight);
      _hist_nch200eta15->fill(numP15, weight);
      _hist_nch200eta30->fill(numP30, weight);
      _hist_nch200eta50->fill(numP50, weight);
      _hist_mean_nch_200->fill(_hist_mean_nch_200->binMean(0), numP50);
    }
    else if (fuzzyEquals(sqrtS, 900.0, 1E-4)) {
      _hist_nch900->fill(numP50, weight);
      _hist_nch900eta05->fill(numP05, weight);
      _hist_nch900eta15->fill(numP15, weight);
      _hist_nch900eta30->fill(numP30, weight);
      _hist_nch900eta50->fill(numP50, weight);
      _hist_mean_nch_900->fill(_hist_mean_nch_900->binMean(0), numP50);
    }
  }
  


  void UA5_1989_S1926373::finalize() {
    // Normalise to area of refhistos
    /// @todo Use generator cross-sections
    normalize(_hist_nch200, 2.011);
    normalize(_hist_nch900, 2.0434);
    normalize(_hist_nch200eta0point5, 1.01255);
    normalize(_hist_nch200eta1point5, 1.0191);
    normalize(_hist_nch200eta3, 1.02615);
    normalize(_hist_nch200eta5, 1.03475);
    normalize(_hist_nch900eta0point5, 1.0035);
    normalize(_hist_nch900eta1point5, 1.01405);
    normalize(_hist_nch900eta3, 1.03055);
    normalize(_hist_nch900eta5, 1.02791);
    // Scale to total number of weights, possible only if binwidth == 1
    scale(_hist_mean_nch_200, 1./sumOfWeights());
    scale(_hist_mean_nch_900, 1./sumOfWeights());
  }



}
