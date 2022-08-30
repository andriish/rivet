// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Xi_c+ spectrum
  class CLEOII_1996_I404590 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOII_1996_I404590);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::abspid==4232), "UFS");
      // book histos
      book(_h_Xi_cPlus,2,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double mH2 =sqr(2.46767);
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
        const double mom = p.momentum().p();
	const double xp = mom/sqrt(0.25*sqr(sqrtS()) - mH2);
	_h_Xi_cPlus->fill(xp);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Xi_cPlus,1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Xi_cPlus;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOII_1996_I404590);

}
