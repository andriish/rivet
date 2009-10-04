// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/TriggerCDFRun0Run1.hh"

namespace Rivet {


  class CDF_1988_S1865951 : public Analysis {
  public:

    /// Constructor
    CDF_1988_S1865951() 
      : Analysis("CDF_1988_S1865951") 
    {
      setBeams(PROTON, ANTIPROTON);
    }


    /// @name Analysis methods
    //@{
    
    /// Book histograms
    void init() {
      addProjection(TriggerCDFRun0Run1(), "Trigger");
      const ChargedFinalState cfs(-1.0, 1.0, 0.4*GeV);
      addProjection(cfs, "CFS");
      addProjection(TotalVisibleMomentum(cfs), "Mom");
      addProjection(Beam(), "Beam");

      _hist_pt1800 = bookHistogram1D(1, 1, 1);
      _hist_pt630 = bookHistogram1D(2, 1, 1);
    }
    
    
    /// Do the analysis
    void analyze(const Event& event) {
      // Trigger
      const bool trigger = applyProjection<TriggerCDFRun0Run1>(event, "Trigger").minBiasDecision();
      if (!trigger) vetoEvent;
      const double weight = event.weight();

      const double sqrtS = applyProjection<Beam>(event, "Beam").sqrtS();
      const FinalState& trackfs = applyProjection<ChargedFinalState>(event, "CFS");
      
      foreach (Particle p, trackfs.particles()) {
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
