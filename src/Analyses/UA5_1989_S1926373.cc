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

    const ChargedFinalState cfs;
    const ChargedFinalState cfs05(-0.5,0.5);
    const ChargedFinalState cfs15(-1.5,1.5);
    const ChargedFinalState cfs30(-3.,3.);
    const ChargedFinalState cfs50(-5.0,5.0);
      
    addProjection(cfs,   "CFSAll");
    addProjection(cfs05, "CFS05");
    addProjection(cfs15, "CFS15");
    addProjection(cfs30, "CFS30");
    addProjection(cfs50, "CFS50");
  }


  // Book histograms
  void UA5_1989_S1926373::init() {
      // _hist_nch200 and _hist_nch200eta5 use the same data but different binning
      // the same is true for _hist_nch900 and _hist_nch900eta5
     _hist_nch200           = bookHistogram1D(1,1,1); 
     _hist_nch900           = bookHistogram1D(2,1,1);
     _hist_nch200eta0point5 = bookHistogram1D(3,1,1);
     _hist_nch200eta1point5 = bookHistogram1D(4,1,1);
     _hist_nch200eta3       = bookHistogram1D(5,1,1);
     _hist_nch200eta5       = bookHistogram1D(6,1,1);
     _hist_nch900eta0point5 = bookHistogram1D(7,1,1);
     _hist_nch900eta1point5 = bookHistogram1D(8,1,1);
     _hist_nch900eta3       = bookHistogram1D(9,1,1);
     _hist_nch900eta5       = bookHistogram1D(10,1,1);
     _hist_mean_nch_200     = bookHistogram1D(11,1,1); 
     _hist_mean_nch_900     = bookHistogram1D(12,1,1);
  } 


  // Do the analysis
    void UA5_1989_S1926373::analyze(const Event& event) {
        Log log = getLog();

        const double sqrtS = applyProjection<Beam>(event, "Beams").sqrtS();
        const double weight = event.weight();

        // Minimum Bias trigger requirements from the hodoscopes
        int n_trig_1 = 0;
        int n_trig_2 = 0;
        
        const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(event, "CFSAll");
        foreach (const Particle& p, cfs.particles()) {
             double eta = p.momentum().pseudorapidity();
             if ( ( -5.6 < eta ) && ( eta < -2.0 ) ) n_trig_1++;
             else if ( ( 2.0 < eta ) && ( eta < 5.6 ) ) n_trig_2++;
        }
        
        // Require at least one coincidence hit in trigger hodoscopes
        if ( n_trig_1* n_trig_2 < 1. ) vetoEvent; 
        getLog() << Log::DEBUG << "Trigger 1: " << n_trig_1 << " Trigger 2: " << n_trig_2 << endl;
        
        // Declare final states in several eta regions 
        const ChargedFinalState& cfs05 = applyProjection<ChargedFinalState>(event, "CFS05");
        const ChargedFinalState& cfs15 = applyProjection<ChargedFinalState>(event, "CFS15");
        const ChargedFinalState& cfs30 = applyProjection<ChargedFinalState>(event, "CFS30");
        const ChargedFinalState& cfs50 = applyProjection<ChargedFinalState>(event, "CFS50");

        // Count particles in final states
        const size_t numP05 = cfs05.particles().size();
        const size_t numP15 = cfs15.particles().size();
        const size_t numP30 = cfs30.particles().size();
        const size_t numP50 = cfs50.particles().size();
        
        // Fill histograms
        if (fuzzyEquals(sqrtS, 200.0, 1E-4)) {
            _hist_nch200->fill(numP50, weight);
            _hist_nch200eta0point5->fill(numP05, weight);
            _hist_nch200eta1point5->fill(numP15, weight);
            _hist_nch200eta3->fill(numP30, weight);
            _hist_nch200eta5->fill(numP50, weight);
            _hist_mean_nch_200->fill(_hist_mean_nch_200->binMean(0), numP50);
            }
        else if (fuzzyEquals(sqrtS, 900.0, 1E-4)) {
            _hist_nch900->fill(numP50, weight);
            _hist_nch900eta0point5->fill(numP05, weight);
            _hist_nch900eta1point5->fill(numP15, weight);
            _hist_nch900eta3->fill(numP30, weight);
            _hist_nch900eta5->fill(numP50, weight);
            _hist_mean_nch_900->fill(_hist_mean_nch_900->binMean(0), numP50);
            }
        }

  void UA5_1989_S1926373::finalize() {
    // Normalise to area of refhistos
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

