// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief rho0 and K*+/- spectra at 35 GeV 
  class JADE_1984_I203145 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(JADE_1984_I203145);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      _h_rho   = bookHisto1D(2, 1, 1);
      _h_kstar = bookHisto1D(3, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get event weight for histo filling
      const double weight = event.weight();

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      foreach (const Particle& p, apply<UnstableParticles>(event, "UFS").
	       particles(Cuts::abspid==323 or Cuts::pid==113)) {
	double xE = p.E()/meanBeamMom;
	double modp = p.p3().mod();
	double beta = modp/p.E();
	if(p.pdgId()==113) {
	  _h_rho->fill(xE,weight/beta);
	}
	else {
	  _h_kstar->fill(xE,weight/beta);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_rho  , sqr(sqrtS())*crossSection()/nanobarn/sumOfWeights());
      scale(_h_kstar, sqr(sqrtS())*crossSection()/nanobarn/sumOfWeights());

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_rho, _h_kstar;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(JADE_1984_I203145);


}
