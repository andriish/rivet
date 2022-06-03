// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief f_1 and f'_1 spectra
  class DELPHI_2003_I628566 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DELPHI_2003_I628566);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_n_f1     , 1, 1, 1);
      book(_n_f1prime, 1, 1, 2);
      book(_h_f1     , 2, 1, 1);
      book(_h_f1prime, 2, 1, 2);
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

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles(Cuts::pid==20223 ||
					     Cuts::pid==20333 )) {
        double xp = p.p3().mod()/meanBeamMom;
	if(p.pid()==20223) {
	  _n_f1     ->fill(xp);
	  _h_f1     ->fill(xp);
	}
	else {
	  _n_f1prime->fill(xp);
	  _h_f1prime->fill(xp);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_n_f1     , 1./sumOfWeights());
      scale(_n_f1prime, 1./sumOfWeights());
      scale(_h_f1     , 1./sumOfWeights());
      scale(_h_f1prime, 1./sumOfWeights());
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _n_f1,_n_f1prime;
    Histo1DPtr _h_f1,_h_f1prime;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(DELPHI_2003_I628566);

}
