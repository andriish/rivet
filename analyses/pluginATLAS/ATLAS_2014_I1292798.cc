// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief chi_c production at 7 TeV
  class ATLAS_2014_I1292798 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2014_I1292798);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      for(unsigned int ichi=0;ichi<2;++ichi) {
	for(unsigned int ipT=0;ipT<2;++ipT) {
	  for(unsigned int iprompt=0;iprompt<2;++iprompt) {
	    book(_h_pT[ichi][iprompt][ipT],1+ichi+2*iprompt+4*ipT,1,1);
	  }
	  book(_h_pT[ichi][2][ipT],"TMP/chi_"+toString(ichi)+"_"+toString(ipT),refData(1+ichi+4*ipT,1,1));
	}
      }
      for(unsigned int iprompt=0;iprompt<2;++iprompt) {
	for(unsigned int ifeed=0;ifeed<2;++ifeed) {
	  book(_h_Jpsi[iprompt][ifeed],"TMP/jpsi_"+toString(iprompt)+"_"+toString(ifeed),
	       refData(9+iprompt,1,1));
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles(Cuts::pid==443 ||
					     Cuts::pid==20443 ||
					     Cuts::pid==445)) {
	// prompt/non-prompt
	bool nonPrompt = p.fromBottom();
	// J/psi as a reference
	if(p.pid()==443) {
	  if(p.absrap()<0.75) _h_Jpsi[nonPrompt][0]->fill(p.perp());
	}
	else if (p.children().size()==2) {
	  Particle Jpsi;
	  if(p.children()[0].pid()==22 && p.children()[1].pid()==443) {
	    Jpsi=p.children()[1];
	  }
	  else if(p.children()[1].pid()==22 && p.children()[0].pid()==443) {
	    Jpsi=p.children()[0];
	  }
	  else
	    continue;
	  if(Jpsi.absrap()>0.75) continue;
	  unsigned int ichi = p.pid()==20443 ? 0 : 1;
	  _h_Jpsi[nonPrompt][1]->fill(Jpsi.perp());
	  _h_pT[ichi][nonPrompt][0]->fill(Jpsi.perp());
	  _h_pT[ichi][nonPrompt][1]->fill(p   .perp());
	  _h_pT[ichi][2        ][0]->fill(Jpsi.perp());
	  _h_pT[ichi][2        ][1]->fill(p   .perp());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // J/psi -> mu+mu- branching ratio
      double br = 0.05961;
      double fact = br*crossSection()/picobarn/sumOfWeights();
      for(unsigned int iprompt=0;iprompt<3;++iprompt) {
	for(unsigned int ichi=0;ichi<2;++ichi) {
	  for(unsigned int ipT=0;ipT<2;++ipT) {
	    scale(_h_pT[ichi][iprompt][ipT],fact);
	  }
	}
      }
      for(unsigned int iprompt=0;iprompt<2;++iprompt) {
	for(unsigned int ifeed=0;ifeed<2;++ifeed) {
	  scale(_h_Jpsi[iprompt][ifeed],fact);
	}
      }
      // Jpsi feed down prompt
      Scatter2DPtr tmp;
      book(tmp,9,1,1);
      efficiency(_h_Jpsi[0][1],_h_Jpsi[0][0],tmp);
      // ratio prompt chi2/chi1
      book(tmp,10,1,1);
      divide(_h_pT[1][0][0],_h_pT[0][0][0],tmp);
      // ratio non-prompt chi2/chi1
      book(tmp,11,1,1);
      divide(_h_pT[1][1][0],_h_pT[0][1][0],tmp);
      // ratio non-prmopt chi1
      book(tmp,12,1,1);
      divide(_h_pT[0][1][1],_h_pT[0][2][1],tmp);
      // ratio non-prmopt chi2
      book(tmp,13,1,1);
      divide(_h_pT[1][1][1],_h_pT[1][2][1],tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pT[2][3][2],_h_Jpsi[2][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ATLAS_2014_I1292798);

}
