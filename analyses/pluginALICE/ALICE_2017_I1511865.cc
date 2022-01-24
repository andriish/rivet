// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/psi, and psi(2s) production at 5.02 and 13 TeV
  class ALICE_2017_I1511865 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2017_I1511865);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      if (isCompatibleWithSqrtS(13000)) {
	book(_h_JPsi_pT,1,1,1);
	book(_h_JPsi_y,2,1,1);
	book(_h_Psi2S_pT,3,1,1);
	book(_h_Psi2S_y ,4,1,1);
	book(_h_JPsi_pT2,"TMP/JPsi_pY",refData(5,1,1));
	book(_h_JPsi_y2 ,"TMP/JPsi_y" ,refData(6,1,1));
      }
      else if(isCompatibleWithSqrtS(5020) ) {
      	book(_h_JPsi_pT,7,1,1);
      	book(_h_JPsi_y,8,1,1);
      }
      else
        throw UserError("Centre-of-mass energy of the given input is neither 5020 nor 13000 GeV.");

	
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over J/Psi
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==443 or
										Cuts::pid==100443)) {
	// rapidity cut
        double absrap = p.absrap();
	if(absrap<2.5 || absrap>4) continue;
	double xp = p.perp();
	// J/Psi
	if(p.pid()==443) {
	  if(_h_JPsi_pT2) {
	    if(xp>30.) continue;
	    _h_JPsi_pT->fill(xp);
	    _h_JPsi_y->fill(absrap);
	    if(xp<=16.) {
	    _h_JPsi_pT2->fill(xp);
	    _h_JPsi_y2->fill(absrap);
	    }
	  }
	  else {
	    if(xp>12.) continue;
	    _h_JPsi_pT->fill(xp);
	    _h_JPsi_y->fill(absrap);
	  }
	}
	// psi(2S)
	else if(_h_Psi2S_pT) {
	  if(xp>16.) continue;
	  _h_Psi2S_pT->fill(xp);
	  _h_Psi2S_y->fill(absrap);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // factor 1/2 due folding +/- rap
      double fact = 0.5*crossSection()/nanobarn/sumOfWeights();
      // factor 1.5 for rapidity range 2.5-4
      scale(_h_JPsi_pT,fact/1.5);
      scale(_h_JPsi_y ,fact);
      if(_h_Psi2S_pT) {
	scale(_h_Psi2S_pT,fact/1.5);
	scale(_h_Psi2S_y ,fact);
      }
      if(_h_JPsi_pT2) {
	scale(_h_JPsi_pT2,fact/1.5);
	Scatter2DPtr tmp;
	book(tmp,5,1,1);
	divide(_h_Psi2S_pT,_h_JPsi_pT2,tmp);
      }
      if(_h_JPsi_y2) {
	scale(_h_JPsi_y2 ,fact);
	Scatter2DPtr tmp;
	book(tmp,6,1,1);
	divide(_h_Psi2S_y,_h_JPsi_y2,tmp);
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_JPsi_pT,_h_JPsi_y,_h_Psi2S_pT,_h_Psi2S_y;
    Histo1DPtr _h_JPsi_pT2,_h_JPsi_y2;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(ALICE_2017_I1511865);

}
