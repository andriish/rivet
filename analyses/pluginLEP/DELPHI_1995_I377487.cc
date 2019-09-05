// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class DELPHI_1995_I377487 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(DELPHI_1995_I377487);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      _h_K0_x  = bookHisto1D( 8, 1, 1);
      _h_K0_xi = bookHisto1D( 9, 1, 1);
      _h_Ks_x  = bookHisto1D(10, 1, 1);

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

      foreach (const Particle& p, ufs.particles(Cuts::pid==130 or
						Cuts::pid==310 or
						Cuts::abspid==323)) {
        double xp = p.p3().mod()/meanBeamMom;
	if(abs(p.pdgId())==323)
	  _h_Ks_x->fill(xp,weight);
	else {
	  _h_K0_x->fill(xp,weight);
	  _h_K0_xi->fill(-log(xp),weight);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_K0_x, 1./sumOfWeights());
      scale(_h_K0_xi, 1./sumOfWeights());
      scale(_h_Ks_x, 1./sumOfWeights());
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_K0_x,_h_K0_xi,_h_Ks_x;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(DELPHI_1995_I377487);


}
