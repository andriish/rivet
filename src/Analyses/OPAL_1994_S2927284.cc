// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "LWH/AIManagedObject.h"
using namespace AIDA;

namespace Rivet {


  /// @brief OPAL charged particle fragmentation functions
  /// @author Peter Richardson
  class OPAL_1994_S2927284 : public Analysis {
  public:

    /// Constructor
    OPAL_1994_S2927284() : Analysis("OPAL_1994_S2927284")
    {}


    /// @name Analysis methods
    //@{

    void analyze(const Event& e) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = applyProjection<FinalState>(e, "FS");

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (fs.particles().size() < 2) {
        MSG_DEBUG("Failed ncharged cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed ncharged cut");
      // Get event weight for histo filling
      const double weight = e.weight();
      
      // Get beams and average beam momentum
      const ParticlePair& beams = applyProjection<Beam>(e, "Beams").beams();
      const double meanBeamMom = ( beams.first.momentum().vector3().mod() +
                                   beams.second.momentum().vector3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      foreach (const Particle& p, fs.particles()) {
	int id = abs(p.pdgId());
	// charged pions
	if(id==PIPLUS) {
	  _histXpPiPlus->fill(p.momentum().vector3().mod(), weight);
	}
	else if(id==KPLUS) {
	  _histXpKPlus->fill(p.momentum().vector3().mod(), weight);
	}
	else if(id==PROTON) {
	  _histXpProton->fill(p.momentum().vector3().mod(), weight);
	}
      }
    }


    void init() {
      // Projections
      addProjection(Beam(), "Beams");
      addProjection(ChargedFinalState(), "FS");

      _histXpPiPlus = bookHistogram1D( 1, 1, 1);
      _histXpKPlus  = bookHistogram1D( 2, 1, 1);
      _histXpProton = bookHistogram1D( 3, 1, 1);
    }


    /// Finalize
    void finalize() {
      scale(_histXpPiPlus,1./sumOfWeights());
      scale(_histXpKPlus ,1./sumOfWeights());
      scale(_histXpProton,1./sumOfWeights());
    }

    //@}

  private:

    AIDA::IHistogram1D *_histXpPiPlus;
    AIDA::IHistogram1D *_histXpKPlus;
    AIDA::IHistogram1D *_histXpProton;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(OPAL_1994_S2927284);

}
