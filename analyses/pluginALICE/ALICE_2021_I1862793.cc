// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Xi_c0,+ at 13 TeV
  class ALICE_2021_I1862793 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2021_I1862793);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(), "UFS");
      // pT distributions
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_pT_Xi[ix],1+ix,1,1);
	book(_h_pT_D[ix],"TMP/h_pT_D_"+toString(ix),refData(3+ix,1,1));
      }
      book(_h_pT_Lambda,"TMP/h_Lambda",refData(5,1,1));
      book(_h_pT_Xi0   ,"TMP/h_Xi0"   ,refData(5,1,1));
      book(_h_pT_Xi_All,"TMP/h_Xi"    ,refData(6,1,1));
      book(_h_pT_Sigma ,"TMP/h_Sigma" ,refData(6,1,1));
      book(_h_sigma[0],7,1,1);
      book(_h_sigma[1],8,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for( const Particle & p : ufs.particles(Cuts::pid==4122 || Cuts::pid==4132 || Cuts::pid==4232 ||
					      Cuts::pid==4222 || Cuts::pid==4212 || Cuts::pid==4112 ||
					      Cuts::pid==421)) {
	// prompt
	if(p.fromBottom()) continue;
	// no mixing
	if(p.children().size()==1 || p.children()[0].abspid()==p.abspid()) continue;
	// rapidity cut
	if(p.absrap()>0.5) continue;
	double pT = p.perp();
	// fill histos for different particles
	// Xi_c0
	if(p.pid()==4132) {
	  _h_pT_Xi[0]->fill(pT);
	  _h_pT_Xi0->fill(pT);
	  _h_pT_Xi_All->fill(pT);
	  if(pT>1.&&pT<12.) _h_sigma[0]->fill(1.);
	  _h_sigma[0]->fill(2.);
	}
	// Xi_c+
	else if(p.pid()==4232) {
	  _h_pT_Xi[1]->fill(pT);
	  _h_pT_Xi_All->fill(pT);
	  _h_sigma[1]->fill(pT);
	}
	// D0 for reference
	else if(p.pid()==421) {
	  _h_pT_D[0]->fill(pT);
	  _h_pT_D[1]->fill(pT);
	}
	// Lambda_c+for reference
	else if(p.pid()==4122) {
	  _h_pT_Lambda->fill(pT);
	}
	// Sigma_c
	else {
	  _h_pT_Sigma->fill(pT);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double factor = crossSection()/microbarn/sumOfWeights();
      for(unsigned int ix=0;ix<2;++ix) {
	scale(_h_pT_Xi[ix],factor);
	scale(_h_pT_D[ix] ,factor);
	// ratio 
	Scatter2DPtr tmp;
	book(tmp,3+ix,1,1);
	divide(_h_pT_Xi[ix],_h_pT_D[ix],tmp);
	scale(_h_sigma[ix],factor);
      }
      // xi_c0/ lambda_c
      Scatter2DPtr tmp;
      book(tmp,5,1,1);
      divide(_h_pT_Xi0,_h_pT_Lambda,tmp);
      // xi_c/sigma_c
      book(tmp,6,1,1);
      divide(_h_pT_Xi_All,_h_pT_Sigma,tmp);
      // remove bin width for xi_c+ total
      scale(_h_sigma[1],8.);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pT_Xi[2],_h_pT_D[2],_h_pT_Lambda,_h_pT_Xi0,_h_pT_Xi_All,_h_pT_Sigma;
    Histo1DPtr _h_sigma[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ALICE_2021_I1862793);

}
