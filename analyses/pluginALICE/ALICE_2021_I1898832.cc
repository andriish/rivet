// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/psi production at 13 TeV (central rapidity)
  class ALICE_2021_I1898832 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2021_I1898832);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_h_pT,1,1,1);
      book(_h_y1,2,1,1);
      book(_h_y2,3,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over J/Psi
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==443)) {
        double absrap = p.absrap();
	if(absrap>0.9) continue;
        double xp = p.perp();
	_h_pT->fill(xp);
	_h_y1->fill(absrap);
	_h_y2->fill(absrap);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/microbarn/sumOfWeights();
      // 1.8 for the rapidity -0.9 to 0.9 in the double differential
      scale(_h_pT,fact/1.8);
      scale(_h_y1,fact);
      // 0.5 as fold + and - rapidity
      scale(_h_y2,0.5*fact);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_pT,_h_y1,_h_y2;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(ALICE_2021_I1898832);

}
