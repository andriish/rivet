// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief eta_c production at 7,8 TeV
  class LHCB_2015_I1316329 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2015_I1316329);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      unsigned int iloc=0;
      if (isCompatibleWithSqrtS(7000)) {
	iloc = 1;
      }
      else if  (isCompatibleWithSqrtS(8000)) {
	iloc = 2;
      }
      else
	throw UserError("Centre-of-mass energy of the given input is neither 7 or 8 TeV.");
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_etac[ix],iloc+2*ix,1,1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles(Cuts::pid==441)) {
	// prompt/non-prompt
	bool nonPrompt = p.fromBottom();
        double absrap = p.absrap();
        double xp = p.perp();
	if(absrap>2. && absrap<4.5) {
	  _h_etac[nonPrompt]->fill(xp);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 1/2 due rapidity folding +/-
      double factor = 0.5*crossSection()/nanobarn/sumOfWeights();
      for(unsigned int ix=0;ix<2;++ix)
	scale(_h_etac[ix],factor);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_etac[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2015_I1316329);

}
