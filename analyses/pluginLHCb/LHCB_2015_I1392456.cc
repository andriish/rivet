// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief Upsilon production at 7,8 TeV
  class LHCB_2015_I1392456 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2015_I1392456);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      int iloc=-1;
      if (isCompatibleWithSqrtS(7000)) {
	iloc = 0;
      }
      else if  (isCompatibleWithSqrtS(8000)) {
	iloc = 10;
      }
      else
	throw UserError("Centre-of-mass energy of the given input is neither 7 or 8 TeV.");

      vector<double> ybins={2.0,2.5,3.0,3.5,4.0,4.5};
      for(unsigned int iups=0;iups<3;++iups) {
	for(unsigned int iy=0;iy<5;++iy) {
	  Histo1DPtr tmp;
	  _h_Ups[iups].add(ybins[iy],ybins[iy+1],book(tmp,1+iups+iloc,1,1+iy));
	}
	book(_h_Ups_pT[iups],4+iloc,1,1+iups);
	book(_h_Ups_y [iups],5+iloc,1,1+iups);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==553 or Cuts::pid==100553 or Cuts::pid==200553)) {
        double absrap = p.absrap();
        double xp = p.perp();
	if(absrap<2. || absrap>4.5) continue;
	unsigned int iups=p.pid()/100000;
	_h_Ups_pT[iups]  ->fill(xp);
	if(xp<30.) _h_Ups_y[iups]   ->fill(absrap);
	_h_Ups[iups]      .fill(absrap,xp);	
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 1/2 due rapidity folding +/-
      double factor = 0.5*crossSection()/picobarn/sumOfWeights();
      // branching ratios for upsilon
      vector<double> br = {0.0248,0.0193,0.0218};
      for(unsigned int iups=0;iups<3;++iups) {
	scale(_h_Ups_pT  [iups],factor*br[iups]);
	scale(_h_Ups_y   [iups],factor/1000.*br[iups]);
	_h_Ups[iups].scale(factor*br[iups],this);
      }
      unsigned int iloc = isCompatibleWithSqrtS(8000) ? 10 : 0;
      for(unsigned int ix=0;ix<_h_Ups[0].histos().size();++ix) {
	Scatter2DPtr tmp;
	for(unsigned int iups=1;iups<3;++iups) {
	  book(tmp,6+iups+iloc,1,1+ix);
	  divide(_h_Ups[iups].histos()[ix],_h_Ups[0].histos()[ix],tmp);
	}
      }
      for(unsigned int iups=1;iups<3;++iups) {
	Scatter2DPtr tmp;
	book(tmp,9+iloc,1,iups);
	divide(_h_Ups_pT[iups],_h_Ups_pT[0],tmp);
	book(tmp,10+iloc,1,iups);
	divide(_h_Ups_y[iups],_h_Ups_y[0],tmp);
      }
      Scatter2DPtr tmp;
      book(tmp,9+iloc,1,3);
      divide(_h_Ups_pT[2],_h_Ups_pT[1],tmp);
      book(tmp,10+iloc,1,3);
      divide(_h_Ups_y[2],_h_Ups_y[1],tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Ups_pT[3], _h_Ups_y[3];
    BinnedHistogram _h_Ups[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2015_I1392456);

}
