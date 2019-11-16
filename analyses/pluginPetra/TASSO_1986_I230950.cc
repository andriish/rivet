// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class TASSO_1986_I230950 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TASSO_1986_I230950);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      _h_p = bookHisto1D(1, 1, 1);
      _h_x = bookHisto1D(2, 1, 1);

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
	       particles(Cuts::pid==PID::PI0)) {
	double modp = p.p3().mod();
	double xE = p.E()/meanBeamMom;
	double beta = modp/p.E();
	_h_p->fill(modp,weight);
	_h_x->fill(xE  ,weight/beta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_p, crossSection()/nanobarn/sumOfWeights());
      scale(_h_x, sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_p, _h_x;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TASSO_1986_I230950);


}
