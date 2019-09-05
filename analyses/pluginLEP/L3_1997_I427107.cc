// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief eta' and omega production at mz
  class L3_1997_I427107 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(L3_1997_I427107);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // // Book histograms
      _histXpOmega   = bookHisto1D( 5, 1, 1);
      _histLnXpOmega = bookHisto1D( 6, 1, 1);
      _histXpEtaP1   = bookHisto1D( 7, 1, 1);
      _histLnXpEtaP1 = bookHisto1D( 8, 1, 1);
      _histXpEtaP2   = bookHisto1D( 9, 1, 1);
      _histLnXpEtaP2 = bookHisto1D(10, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.                                                    
      const FinalState& fs = apply<FinalState>(event, "FS");
      if (fs.particles().size() < 2) {
      	MSG_DEBUG("Failed ncharged cut");
      	vetoEvent;
      }
      MSG_DEBUG("Passed ncharged cut");

      // Get event weight for histo filling                                                                                                    
      const double weight = event.weight();

      // Get beams and average beam momentum                                                                                                
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      // Final state of unstable particles to get particle spectra
      const Particles& mesons = apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==PID::ETAPRIME or
      										 Cuts::pid==PID::OMEGA);

      foreach (const Particle& p, mesons) {
      	double xp = p.p3().mod()/meanBeamMom;
	double xi = log(1./xp);
      	if(p.pdgId()==PID::ETAPRIME) {
      	  _histXpEtaP1->fill(xp, weight);
      	  _histLnXpEtaP1->fill(xi, weight);
      	  _histXpEtaP2->fill(xp, weight);
      	  _histLnXpEtaP2->fill(xi, weight);
      	}
      	else {
      	  _histXpOmega->fill(xp, weight);
      	  _histLnXpOmega->fill(xi, weight);
      	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_histXpEtaP1  , 1./sumOfWeights());
      scale(_histLnXpEtaP1, 1./sumOfWeights());
      scale(_histXpEtaP2  , 1./sumOfWeights());
      scale(_histLnXpEtaP2, 1./sumOfWeights());
      scale(_histXpOmega  , 1./sumOfWeights());
      scale(_histLnXpOmega, 1./sumOfWeights());
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _histXpEtaP1,_histXpEtaP2;
    Histo1DPtr _histLnXpEtaP1,_histLnXpEtaP2;
    Histo1DPtr _histXpOmega;
    Histo1DPtr _histLnXpOmega;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(L3_1997_I427107);


}
