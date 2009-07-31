#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analyses/UA5_1982_S875503.hh"

namespace Rivet {
    
    /// Default constructor
    UA5_1982_S875503::UA5_1982_S875503()
   : Analysis("UA5_1982_S875503") {
      const FinalState fs;
      addProjection(Beam(), "Beam");
      addProjection(fs, "FS");
    }
 
    /// @name Analysis methods
    //@{
    void UA5_1982_S875503::init() 
    { 
    _hist_etapp = bookHistogram1D(3,1,1);
    _hist_etappbar = bookHistogram1D(4,1,1);
    }
    
    void UA5_1982_S875503::analyze(const Event& event) {
      const Beam b = applyProjection<Beam>(event, "Beam");
      const FinalState& fs = applyProjection<FinalState>(event, "FS");
      const double weight = event.weight();
      foreach (const Particle& p, fs.particles())
      {
        if ( b.beams().first.pdgId() == b.beams().second.pdgId())
        {
          _hist_etapp->fill(p.momentum().pseudorapidity(), weight);
        }
        else if ( b.beams().first.pdgId() != b.beams().second.pdgId())
        {
          _hist_etappbar->fill(p.momentum().pseudorapidity(), weight);
        }
      }

    }
    
    void UA5_1982_S875503::finalize() {
      normalize(_hist_etapp, 5.28);
      normalize(_hist_etappbar, 5.29);
    }
    //@}

  }

