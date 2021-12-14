// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {
  class CMS_2015_I1342266 : public Analysis {
  public:
    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2015_I1342266);

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      {Histo1DPtr tmp; _h_Y1s_dy_pT.add(0.0, 0.6, book(tmp, 1, 1, 1));}
      {Histo1DPtr tmp; _h_Y1s_dy_pT.add(0.6, 1.2, book(tmp, 2, 1, 1));}
      book(_h_Y1s_pT, 3, 1, 1);
      
      {Histo1DPtr tmp; _h_Y2s_dy_pT.add(0.0, 0.6, book(tmp, 1, 1, 2));}
      {Histo1DPtr tmp; _h_Y2s_dy_pT.add(0.6, 1.2, book(tmp, 2, 1, 2));}
      book(_h_Y2s_pT, 3, 1, 2);

      {Histo1DPtr tmp; _h_Y3s_dy_pT.add(0.0, 0.6, book(tmp, 1, 1, 3));}
      {Histo1DPtr tmp; _h_Y3s_dy_pT.add(0.6, 1.2, book(tmp, 2, 1, 3));}
      book(_h_Y3s_pT, 3, 1, 3);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles(Cuts::pid==553 || Cuts::pid==100553 || Cuts::pid==200553)) {
        double absrap = p.absrap();
        double xp = p.perp();
	if(absrap>1.2) continue;
        if(p.pid() == 553) {
	  _h_Y1s_pT  ->fill(xp);
	  _h_Y1s_dy_pT.fill(absrap,xp);
	}
        else if (p.pid() == 100553) {
	  _h_Y2s_pT  ->fill(xp);
	  _h_Y2s_dy_pT.fill(absrap,xp);
	}
        else if (p.pid() == 200553) {
	  _h_Y3s_pT  ->fill(xp);
	  _h_Y3s_dy_pT.fill(absrap,xp);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double factor = crossSection() / femtobarn/ sumOfWeights();
      _h_Y1s_dy_pT.scale(factor * 0.0248 , this); //branching ratio
      _h_Y2s_dy_pT.scale(factor * 0.0193 , this);
      _h_Y3s_dy_pT.scale(factor * 0.0218 , this);

      scale(_h_Y1s_pT,factor * 0.0248);
      scale(_h_Y2s_pT,factor * 0.0193);
      scale(_h_Y3s_pT,factor * 0.0218);
      for(unsigned int i=1;i<3;++i) {
	Scatter2DPtr tmp;
	book(tmp,i+3,1,1);
	divide(_h_Y2s_dy_pT.histos()[i-1],_h_Y1s_dy_pT.histos()[0],tmp);
	book(tmp,i+3,1,2);
	divide(_h_Y3s_dy_pT.histos()[i-1],_h_Y1s_dy_pT.histos()[0],tmp);
      }
      Scatter2DPtr tmp;
      book(tmp,6,1,1);
      divide(_h_Y2s_pT,_h_Y1s_pT,tmp);
      book(tmp,6,1,2);
      divide(_h_Y3s_pT,_h_Y1s_pT,tmp);
    }
    
  private:
    BinnedHistogram _h_Y1s_dy_pT;
    BinnedHistogram _h_Y2s_dy_pT;
    BinnedHistogram _h_Y3s_dy_pT;
    Histo1DPtr _h_Y1s_pT;
    Histo1DPtr _h_Y2s_pT;
    Histo1DPtr _h_Y3s_pT;
  };
  DECLARE_RIVET_PLUGIN(CMS_2015_I1342266);
}
