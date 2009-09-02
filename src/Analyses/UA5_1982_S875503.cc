// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {

  class UA5_1982_S875503 : public Analysis {
  public:
    
    /// Default constructor
    UA5_1982_S875503()
      : Analysis("UA5_1982_S875503") 
    {
      const ChargedFinalState cfs(-3.5, 3.5);
      addProjection(Beam(), "Beam");
      addProjection(cfs, "CFS");
    }
  

    /// @name Analysis methods
    //@{

    void init() 
    { 
      _hist_nch_pp    = bookHistogram1D(2,1,1);
      _hist_nch_ppbar = bookHistogram1D(2,1,2);
      _hist_eta_pp    = bookHistogram1D(3,1,1);
      _hist_eta_ppbar = bookHistogram1D(4,1,1);
    }
    
    
    void analyze(const Event& event) {
      const Beam b = applyProjection<Beam>(event, "Beam");
      const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(event, "CFS");
      const double weight = event.weight();
      
      // Different trigger implementations for ppbar and pp!
      int n_trig_1(0), n_trig_2(0);
      foreach (const Particle& p, cfs.particles()) {
        double eta = p.momentum().pseudorapidity();
        if (inRange(eta, -5.6, -2.0)) n_trig_1 += 1;
        else if (inRange(eta, 2.0, 5.6)) n_trig_2 += 1;
      }
      
      // Trigger requirements
      const bool samebeam = (b.beams().first.pdgId() == b.beams().second.pdgId());
      if (samebeam) {
        // PP
        if (n_trig_1 == 0 || n_trig_2 == 0) vetoEvent; 
      } else {
        // PPbar
        /// @todo Is this actually the exact trigger requirement?
        if (n_trig_1 * n_trig_2 < 4) vetoEvent;
      }
      
      // Iterate over all FS particles and fill histograms
      foreach (const Particle& p, cfs.particles()) {
        if (samebeam) {
          // PP collision
          _hist_eta_pp->fill(fabs(p.momentum().pseudorapidity()), weight);
        } else {
          // PPbar collision
          _hist_eta_ppbar->fill(fabs(p.momentum().pseudorapidity()), weight);
        }
      }
      
      // Fill mean charged multiplicity histos
      if (samebeam) {
        // PP
        _hist_nch_pp->fill(_hist_nch_pp->binMean(0), cfs.particles().size());
      } else {
        // PPbar 
        _hist_nch_ppbar->fill(_hist_nch_ppbar->binMean(0), cfs.particles().size());
      }
      
    }
    
    
    void finalize() {
      scale(_hist_nch_pp,    1./sumOfWeights());
      scale(_hist_nch_ppbar, 1./sumOfWeights());
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
