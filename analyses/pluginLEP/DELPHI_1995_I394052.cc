// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief K+/- and protno spectra
  class DELPHI_1995_I394052 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(DELPHI_1995_I394052);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");

      // Book histograms
      _h_kaon_p    = bookHisto1D(3, 1, 1);
      _h_kaon_x    = bookHisto1D(5, 1, 1);
      _h_proton_p  = bookHisto1D(4, 1, 1);
      _h_proton_x  = bookHisto1D(6, 1, 1);
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
      foreach (const Particle& p, fs.particles(Cuts::abspid==321 or
					       Cuts::abspid==2212)) {
	double modp = p.p3().mod();
        double xp = modp/meanBeamMom;
	if(abs(p.pdgId())==321) {
	  _h_kaon_p  ->fill(modp,weight);
	  _h_kaon_x  ->fill(xp  ,weight);
	}
	else {
	  _h_proton_p->fill(modp,weight);
	  _h_proton_x->fill(xp  ,weight);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_kaon_p  ,1./sumOfWeights());
      scale(_h_kaon_x  ,1./sumOfWeights());
      scale(_h_proton_p,1./sumOfWeights());
      scale(_h_proton_x,1./sumOfWeights());
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_kaon_p, _h_kaon_x, _h_proton_p, _h_proton_x;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(DELPHI_1995_I394052);


}
