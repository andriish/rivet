// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/TriggerUA5.hh"

namespace Rivet {

  class UA5_1982_S875503 : public Analysis {
  public:
    
    /// Default constructor
    UA5_1982_S875503() : Analysis("UA5_1982_S875503") {
      //
    }
  

    /// @name Analysis methods
    //@{

    void init() { 
      addProjection(TriggerUA5(), "Trigger");
      addProjection(ChargedFinalState(-3.5, 3.5), "CFS");

      _hist_nch_pp    = bookHistogram1D(2,1,1);
      _hist_nch_ppbar = bookHistogram1D(2,1,2);
      _hist_eta_pp    = bookHistogram1D(3,1,1);
      _hist_eta_ppbar = bookHistogram1D(4,1,1);
    }
    
    
    void analyze(const Event& event) {
      // Trigger
      const TriggerUA5& trigger = applyProjection<TriggerUA5>(event, "Trigger");
      if (!trigger.nsdDecision()) vetoEvent;
      const double weight = event.weight(); 
      
      // Iterate over all tracks and fill histograms
      const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(event, "CFS");
      foreach (const Particle& p, cfs.particles()) {
        if (trigger.samebeams()) {
          // PP collision
          _hist_eta_pp->fill(fabs(p.momentum().pseudorapidity()), weight);
        } else {
          // PPbar collision
          _hist_eta_ppbar->fill(fabs(p.momentum().pseudorapidity()), weight);
        }
      }
      
      // Fill mean charged multiplicity histos
      if (trigger.samebeams()) {
        // PP
        _hist_nch_pp->fill(_hist_nch_pp->binMean(0), cfs.particles().size());
      } else {
        // PPbar 
        _hist_nch_ppbar->fill(_hist_nch_ppbar->binMean(0), cfs.particles().size());
      }
      
    }
    
    
    void finalize() {
      scale(_hist_nch_pp,    1.0/sumOfWeights());
      scale(_hist_nch_ppbar, 1.0/sumOfWeights());
      normalize(_hist_eta_pp,    5.28);
      normalize(_hist_eta_ppbar, 5.29);
    }

    //@}
    
  
  private:
    
    /// @name Histogram collections
    //@{
    AIDA::IHistogram1D* _hist_nch_pp;
    AIDA::IHistogram1D* _hist_nch_ppbar;
    AIDA::IHistogram1D* _hist_eta_pp;
    AIDA::IHistogram1D* _hist_eta_ppbar;
    //@}
    
  };
  
  
  
  // This global object acts as a hook for the plugin system
  AnalysisBuilder<UA5_1982_S875503> plugin_UA5_1982_S875503;
  
}
