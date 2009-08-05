#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/Analyses/UA5_1982_S875503.hh"

namespace Rivet {
    
    /// Default constructor
    UA5_1982_S875503::UA5_1982_S875503()
   : Analysis("UA5_1982_S875503") {
      const ChargedFinalState cfs(-3.5, 3.5);
      addProjection(Beam(), "Beam");
      addProjection(cfs, "CFS");
    }
 
    /// @name Analysis methods
    //@{
    void UA5_1982_S875503::init() 
    { 
    _hist_nch_pp    = bookHistogram1D(2,1,1);
    _hist_nch_ppbar = bookHistogram1D(2,1,2);
    _hist_eta_pp    = bookHistogram1D(3,1,1);
    _hist_eta_ppbar = bookHistogram1D(4,1,1);
    }
    
    void UA5_1982_S875503::analyze(const Event& event) {
      const Beam b = applyProjection<Beam>(event, "Beam");
      const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(event, "CFS");
      const double weight = event.weight();
      int b_pdg = b.beams().first.pdgId() * b.beams().second.pdgId();

      // Different trigger implementations for ppbar and pp!
      int n_trig_1 = 0;
      int n_trig_2 = 0;

      foreach (const Particle& p, cfs.particles()) {
           double eta = p.momentum().pseudorapidity();
           if ( ( -5.6 < eta ) && ( eta < -2.0 ) ) n_trig_1++;
           else if ( ( 2.0 < eta ) && ( eta < 5.6 ) ) n_trig_2++;
      }
      // PP first
      if ( ( b_pdg > 0. ) && ( n_trig_1* n_trig_2 < 1. ) ) {
          vetoEvent; 
      }
      // PPbar trigger requirements
      else if ( ( b_pdg < 0. ) && ( n_trig_1* n_trig_2 < 4. ) ) {
          vetoEvent;
      }


      // Iterate over all FS particles and fill histograms
      foreach (const Particle& p, cfs.particles()) {
        // PP collision
        if ( b.beams().first.pdgId() == b.beams().second.pdgId()) {
            _hist_eta_pp->fill(fabs(p.momentum().pseudorapidity()), weight);
            }
        // PPbar collision
        else if ( b.beams().first.pdgId() != b.beams().second.pdgId()) {
            _hist_eta_ppbar->fill(fabs(p.momentum().pseudorapidity()), weight);
            }
      }

      // Fill mean charged multiplicity histos
      // PP first
      if ( b_pdg > 0. ) {
          _hist_nch_pp->fill(_hist_nch_pp->binMean(0), cfs.particles().size());
      }
      // PPbar 
      else if ( b_pdg < 0. ) {
          _hist_nch_ppbar->fill(_hist_nch_ppbar->binMean(0), cfs.particles().size());
      }

    }
    
    void UA5_1982_S875503::finalize() {
      scale(_hist_nch_pp,    1./sumOfWeights());
      scale(_hist_nch_ppbar, 1./sumOfWeights());
      normalize(_hist_eta_pp,    5.28);
      normalize(_hist_eta_ppbar, 5.29);
    }
    //@}

  }

