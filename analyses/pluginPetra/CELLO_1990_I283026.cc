// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief Strange hadrons at 35GeV
  class CELLO_1990_I283026 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CELLO_1990_I283026);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      _h_K0     = bookHisto1D(1, 1, 1);
      _h_Kstar  = bookHisto1D(2, 1, 1);
      _h_Lambda = bookHisto1D(3, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();
      
      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");
      
      // Get event weight for histo filling
      const double weight = event.weight();
      
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
      
      foreach (const Particle& p, ufs.particles()) {
        const int id = p.abspid();
        if (id == PID::K0S || id == PID::K0L) {
          double xE = p.E()/meanBeamMom;
          _h_K0->fill(xE, weight);
        }
	else if(abs(id)==323) {
          double xE = p.E()/meanBeamMom;
          _h_Kstar->fill(xE, weight);
	}
	else if(abs(id)==3122) {
          double xE = p.E()/meanBeamMom;
          _h_Lambda->fill(xE, weight);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_K0    , 1./sumOfWeights());
      scale(_h_Kstar , 1./sumOfWeights());
      scale(_h_Lambda, 1./sumOfWeights());
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_K0, _h_Kstar, _h_Lambda;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CELLO_1990_I283026);


}
