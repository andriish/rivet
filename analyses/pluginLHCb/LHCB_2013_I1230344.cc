// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief J/psi and Upslion production at 8 TeV
  class LHCB_2013_I1230344 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2013_I1230344);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      vector<double> ybins={2.0,2.5,3.0,3.5,4.0,4.5};
      for(unsigned iy=0;iy<2;++iy) 
	book(_h_Jpsi_y[iy],2,1,iy+1);
      for(unsigned int ix=0;ix<5;++ix) {
	for(unsigned iy=0;iy<2;++iy) {
	  Histo1DPtr tmp;
	  _h_Jpsi[iy].add(ybins[ix],ybins[ix+1],book(tmp,3+ix,1,1+iy));
	}
	Histo1DPtr tmp;
	_h_Jpsi[2].add(ybins[ix],ybins[ix+1],book(tmp,"TMP/Jpsi_"+toString(ix),refData(3+ix,1,1)));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	book(_h_Ups_pT  [iy], 9,1,iy+1);
	book(_h_Ups_y   [iy],10,1,iy+1);
	book(_h_Ups_pT_r[iy],"TMP/Ups_pT_"+toString(iy),refData(16,1,1));
	book(_h_Ups_y_r [iy],"TMP/Ups_y_" +toString(iy),refData(17,1,1));
	for(unsigned int ix=0;ix<5;++ix) {
	  Histo1DPtr tmp;
	  _h_Ups[iy].add(ybins[ix],ybins[ix+1],book(tmp,11+ix,1,1+iy));
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // J/psi
      for (const Particle& p : ufs.particles(Cuts::pid==443)) {
	// prompt/non-prompt
	bool nonPrompt = p.fromBottom();
        double absrap = p.absrap();
        double xp = p.perp();
	if(absrap<2. || absrap>4.5 ||  xp>14.) continue;
	_h_Jpsi_y[nonPrompt]->fill(absrap);
	_h_Jpsi[nonPrompt]   .fill(absrap,xp);
	_h_Jpsi[2]   .fill(absrap,xp);
      }
      // upsilon
      for (const Particle& p : ufs.particles(Cuts::pid==553 or Cuts::pid==100553 or Cuts::pid==200553)) {
        double absrap = p.absrap();
        double xp = p.perp();
	if(absrap<2. || absrap>4.5 ||  xp>15.) continue;
	unsigned int iups=p.pid()/100000;
	_h_Ups_pT[iups]  ->fill(xp);
	_h_Ups_y[iups]   ->fill(absrap);
	_h_Ups_pT_r[iups]->fill(xp);
	_h_Ups_y_r[iups] ->fill(absrap);
	_h_Ups[iups]      .fill(absrap,xp);	
      }
      
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 1/2 due rapidity folding +/-
      double factor = 0.5*crossSection()/nanobarn/sumOfWeights();
      // branching ratios for upsilon
      vector<double> br = {0.0248,0.0193,0.0218};
      for(unsigned int ix=0;ix<2;++ix) {
	scale(_h_Jpsi_y[ix],factor);
      }
      for(unsigned int ix=0;ix<3;++ix) {
	scale(_h_Ups_pT  [ix],factor*br[ix]);
	scale(_h_Ups_y   [ix],factor*br[ix]);
	scale(_h_Ups_pT_r[ix],factor*br[ix]);
	scale(_h_Ups_y_r [ix],factor*br[ix]);
	_h_Ups[ix].scale(factor*br[ix],this);
	_h_Jpsi[ix].scale(factor,this);
      }
      for(unsigned int ix=0;ix<_h_Jpsi[0].histos().size();++ix) {
	Scatter2DPtr tmp;
	book(tmp,3+ix,1,3);
	divide(_h_Jpsi[1].histos()[ix],_h_Jpsi[2].histos()[ix],tmp);
	tmp->scaleY(100.);
      }
      for(unsigned int ix=0;ix<2;++ix) {
	Scatter2DPtr tmp;
	book(tmp,16,1,1+ix);
	divide(_h_Ups_pT_r[ix+1],_h_Ups_pT_r[0],tmp);
	book(tmp,17,1,1+ix);
	divide(_h_Ups_y_r[ix+1],_h_Ups_y_r[0],tmp);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Jpsi_y[2],_h_Ups_pT[3],_h_Ups_y[3],_h_Ups_pT_r[3],_h_Ups_y_r[3];
    BinnedHistogram _h_Jpsi[3],_h_Ups[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2013_I1230344);

}
