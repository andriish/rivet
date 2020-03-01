// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class ARGUS_1993_I357132 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ARGUS_1993_I357132);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      // book histos
      book(_h_x,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const int id2625 = 4124;
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double Emax = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      const double Pmax = sqrt(sqr(Emax)-sqr(2.625));
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==id2625)) {
	double xp = p.momentum().p3().mod()/Pmax;
	_h_x->fill(xp);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_x);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_x;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(ARGUS_1993_I357132);

}
