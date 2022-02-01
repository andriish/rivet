// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Xi_c0 at 7 TeV
  class ALICE_2018_I1642729 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2018_I1642729);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(), "UFS");
      // histos
      book(_h_Xi_pT[0],1,1,1);
      book(_h_Xi_pT[1],"TMP/pT_Xi",refData(2,1,1));
      book(_h_D0_pT   ,"TMP/pT_D0",refData(2,1,1));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for( const Particle & p : ufs.particles(Cuts::pid==4132 || Cuts::pid==421)) {
	// no mixing
	if(p.children().size()==1 || p.children()[0].abspid()==p.abspid()) continue;
	// rapidity cut
	if(p.absrap()>0.5) continue;
	double pT = p.perp();
	if(p.pid()==421) {
	  _h_D0_pT->fill(pT);
	}
	else {
	  _h_Xi_pT[0]->fill(pT);
	  _h_Xi_pT[1]->fill(pT);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double factor = crossSection()/microbarn/sumOfWeights();
      // from PDG value of Xi_x0-> Xi- pi+ and ALICe measurement of ratio to semileptonic
      double br = 0.0197;
      for(unsigned int ix=0;ix<2;++ix) {
	scale(_h_Xi_pT[ix],br*factor);
      }
      scale(_h_D0_pT,factor);
      Scatter2DPtr tmp;
      book(tmp,2,1,1);
      divide(_h_Xi_pT[1],_h_D0_pT,tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Xi_pT[2],_h_D0_pT;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ALICE_2018_I1642729);

}
