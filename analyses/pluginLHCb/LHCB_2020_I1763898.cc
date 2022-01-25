// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief eta_c production at 13 TeV
  class LHCB_2020_I1763898 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2020_I1763898);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_etac[ix],2+2*ix,1,1);
	book(_h_Jpsi[ix],"TMP/Jpsi"+toString(ix),refData(2+2*ix,1,1));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles(Cuts::pid==441|| Cuts::pid==443)) {
	// prompt/non-prompt
	bool nonPrompt = p.fromBottom();
        double absrap = p.absrap();
        double xp = p.perp();
	if(absrap>2. && absrap<4.5) {
	  if(p.pid()==441)
	    _h_etac[nonPrompt]->fill(xp);
	  else
	    _h_Jpsi[nonPrompt]->fill(xp);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 1/2 due rapidity folding +/-
      double factor = 0.5*crossSection()/nanobarn/sumOfWeights();
      for(unsigned int ix=0;ix<2;++ix) {
	scale(_h_etac[ix],factor);
	scale(_h_Jpsi[ix],factor);
	Scatter2DPtr tmp;
	book(tmp,1+2*ix,1,1);
	divide(_h_etac[ix],_h_Jpsi[ix],tmp);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_etac[2],_h_Jpsi[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2020_I1763898);

}
