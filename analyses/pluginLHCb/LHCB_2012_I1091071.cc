// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief Upslion production at 7 TeV
  class LHCB_2012_I1091071 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2012_I1091071);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      vector<double> ybins={2.0,2.5,3.0,3.5,4.0,4.5};
      for(unsigned int ix=0;ix<3;++ix) {
	for(unsigned int iy=0;iy<5;++iy) {
	  Histo1DPtr tmp;
	  _h_Ups[ix].add(ybins[iy],ybins[iy+1],book(tmp,2+iy+5*ix,1,1));
	}
	book(_h_Ups_pT[ix],"TMP/Ups_"+toString(ix),refData(17,1,1));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // upsilon
      for (const Particle& p : ufs.particles(Cuts::pid==553 or Cuts::pid==100553 or Cuts::pid==200553)) {
        double absrap = p.absrap();
        double xp = p.perp();
	if(absrap<2. || absrap>4.5) continue;
	unsigned int iups=p.pid()/100000;
	_h_Ups_pT[iups]  ->fill(xp);
	_h_Ups[iups]      .fill(absrap,xp);	
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 1/2 due rapidity folding +/-
      double factor = 0.5*crossSection()/picobarn/sumOfWeights();
      // branching ratios for upsilon
      vector<double> br = {0.0248,0.0193,0.0218};
      for(unsigned int ix=0;ix<3;++ix) {
	scale(_h_Ups_pT  [ix],factor*br[ix]);
	_h_Ups[ix].scale(factor*br[ix],this);
      }
      for(unsigned int ix=0;ix<2;++ix) {
	Scatter2DPtr tmp;
	book(tmp,17,1,1+ix);
	divide(_h_Ups_pT[ix+1],_h_Ups_pT[0],tmp);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Ups_pT[3];
    BinnedHistogram _h_Ups[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2012_I1091071);

}
