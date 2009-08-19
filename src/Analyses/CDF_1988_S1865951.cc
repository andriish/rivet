#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Analyses/CDF_1988_S1865951.hh"

namespace Rivet {

    /// Default constructor
    CDF_1988_S1865951::CDF_1988_S1865951() 
    : Analysis("CDF_1988_S1865951") {
      const ChargedFinalState cfs(-1.,1., 0.4*GeV);
      addProjection(cfs, "CFS");
      addProjection(ChargedFinalState(-5.9, 5.9), "CFSAll");
      addProjection(TotalVisibleMomentum(cfs), "Mom");
      addProjection(Beam(), "Beam");
    }

    /// @name Analysis methods
    //@{
    void CDF_1988_S1865951::init() { 
      _hist_pt1800 = bookHistogram1D(1, 1, 1);
      _hist_pt630 = bookHistogram1D(2, 1, 1);
    }
    
    void CDF_1988_S1865951::analyze(const Event& event) {
      const double sqrtS = applyProjection<Beam>(event, "Beam").sqrtS();
      const FinalState& fs = applyProjection<ChargedFinalState>(event, "CFS");
      const double weight = event.weight();

        // Minimum Bias trigger requirements from the BBC counters
        int n_trig_1 = 0;
        int n_trig_2 = 0;

        // Event selection based on tracks in VTPC (time projection chambers)
        // Require at least 4 tracks with at least one in each of the forward
        // and backward hemispheres
        int n_backward = 0;
        int n_forward = 0;
        
        const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(event, "CFSAll");
        foreach (const Particle& p, cfs.particles()) {
             double eta = p.momentum().pseudorapidity();
             if (inRange(eta, -5.9, -3.2)) n_trig_1++;
             else if (inRange(eta, 3.2, 5.9)) n_trig_2++;

             if (inRange(eta, -3.0, 0.0)) n_backward++;
             else if (inRange(eta, 0.0, 3.0)) n_forward++;
        }
            
        // Require at least one coincidence hit in both BBC counters
        if (n_trig_1*n_trig_2 < 1) vetoEvent; 
        getLog() << Log::DEBUG << "Trigger 1: " << n_trig_1 << " Trigger 2: " << n_trig_2 << endl;
      
        // Further event selection cut
        if ( (n_backward+n_forward < 4) || (n_backward*n_forward < 1) ) vetoEvent;
        getLog() << Log::DEBUG << " Num. forward: " << n_forward  << " Num. backward: " << n_backward << endl;
      
      foreach (Particle p, fs.particles())
      {
        double pt = p.momentum().pT();
        if (fuzzyEquals(sqrtS, 630/GeV))
        {
          // weight / ( Delta eta * 2pi * pt ) with Delta eta = 2
          _hist_pt630->fill(pt, weight/(2.*2.*M_PI*pt));
        }
        else if (fuzzyEquals(sqrtS, 1800/GeV))
        {
          _hist_pt1800->fill(pt, weight/(2.*2.*M_PI*pt));
        }
      
      }
    }
    
    void CDF_1988_S1865951::finalize() {
    ///@todo Total cross section hard-coded, needs a way to pass variable from pythia.
      scale(_hist_pt630, 32.6/sumOfWeights());
      scale(_hist_pt1800, 38.5/sumOfWeights());
    }
    //@}
    


  }
