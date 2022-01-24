// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  J/psi at 5.02 TeV
  class ALICE_2019_I1735351 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2019_I1735351);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_h_JPsi_pT,2,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over J/Psi
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==443)) {
	// rapidity cut
	if(p.absrap()>0.9) continue;
	double xp = p.perp();
	_h_JPsi_pT->fill(xp);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // factor 1/2 due folding +/- rap and 1.8 due y range
      double fact = 0.5/1.8*crossSection()/nanobarn/sumOfWeights();
      scale(_h_JPsi_pT,fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_JPsi_pT;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ALICE_2019_I1735351);

}
