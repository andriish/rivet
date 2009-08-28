// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  class CDF_1988_S1865951 : public Analysis {

  public:

    /// @name Constructor etc.
    //@{

    /// Constructor
    CDF_1988_S1865951() 
      : Analysis("CDF_1988_S1865951") 
    {
      const ChargedFinalState cfs(-1.,1., 0.4*GeV);
      addProjection(cfs, "CFS");
      addProjection(ChargedFinalState(-5.9, 5.9), "CFSAll");
      addProjection(TotalVisibleMomentum(cfs), "Mom");
      addProjection(Beam(), "Beam");
    }
    
    //@}


    /// @name Analysis methods
    //@{
    
    /// Book histograms
    void init() { 
      _hist_pt1800 = bookHistogram1D(1, 1, 1);
      _hist_pt630 = bookHistogram1D(2, 1, 1);
    }
    
    
    /// Do the analysis
    void analyze(const Event& event) {
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
        if (inRange(eta, -5.9, -3.2)) n_trig_1 += 1;
        else if (inRange(eta, 3.2, 5.9)) n_trig_2 += 1;
        
        if (inRange(eta, -3.0, 0.0)) n_backward += 1;
        else if (inRange(eta, 0.0, 3.0)) n_forward += 1;
      }
      
      // Require at least one coincidence hit in both BBC counters
      if (!n_trig_1 || !n_trig_2) vetoEvent; 
      getLog() << Log::DEBUG << "Trigger 1: " << n_trig_1 << " Trigger 2: " << n_trig_2 << endl;
      
      // Further event selection cut
      if ( (n_backward+n_forward < 4) || (n_backward*n_forward < 1) ) vetoEvent;
      getLog() << Log::DEBUG << " Num. forward: " << n_forward  << " Num. backward: " << n_backward << endl;
      
      foreach (Particle p, fs.particles()) {
        const double pt = p.momentum().pT();
        // Effective weight for d3sig/dp3 = weight / ( Delta eta * 2pi * pt ), with Delta(eta) = 2
        const double eff_weight = weight/(2*TWOPI*pt);
        if (fuzzyEquals(sqrtS, 630/GeV)) {
          _hist_pt630->fill(pt, eff_weight);
        } else if (fuzzyEquals(sqrtS, 1800/GeV)) {
          _hist_pt1800->fill(pt, eff_weight);
        }
        
      }
    }
    
    
    /// Scale histos
    void finalize() {
      /// @todo Total cross section hard-coded, needs a way to pass variable from generator
      scale(_hist_pt630, 32.6/sumOfWeights());
      scale(_hist_pt1800, 38.5/sumOfWeights());
    }
   
    //@}

  private:
    
    /// @name Histos
    //@{
    AIDA::IHistogram1D* _hist_pt630;
    AIDA::IHistogram1D* _hist_pt1800;
    //@}

  };
 
  

  // This global object acts as a hook for the plugin system
  AnalysisBuilder<CDF_1988_S1865951> plugin_CDF_1988_S1865951;

}
