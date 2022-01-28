// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Upsilon production at 2.76 TeV
  class LHCB_2014_I1280929 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2014_I1280929);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      for(unsigned int iups=0;iups<3;++iups) {
	book(_h_Ups_pT [iups],3+iups,1,1);
	book(_h_Ups_y  [iups],6+iups,1,1);
	book(_h_Ups_y_r[iups],"TMP/UPS_y_"+toString(iups),refData(10,1,1));

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
	if(absrap<2. || absrap>4.5 ||  xp>15.) continue;
	unsigned int iups=p.pid()/100000;
	_h_Ups_pT[iups]  ->fill(xp);
	_h_Ups_y[iups]   ->fill(absrap);
	_h_Ups_y_r[iups] ->fill(absrap);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 1/2 due rapidity folding +/-
      double factor = 0.5*crossSection()/nanobarn/sumOfWeights();
      // branching ratios for upsilon
      vector<double> br = {0.0248,0.0193,0.0218};
      for(unsigned int ix=0;ix<3;++ix) {
	scale(_h_Ups_pT  [ix],factor*br[ix]);
	scale(_h_Ups_y   [ix],factor*br[ix]); 
	scale(_h_Ups_y_r [ix],factor*br[ix]);
      }
      for(unsigned int ix=0;ix<2;++ix) {
	Scatter2DPtr tmp;
	book(tmp,9,1,1+ix);
	divide(_h_Ups_pT[ix+1],_h_Ups_pT[0],tmp);
	book(tmp,10,1,1+ix);
	divide(_h_Ups_y_r[ix+1],_h_Ups_y_r[0],tmp);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Ups_pT[3],_h_Ups_y[3],_h_Ups_y_r[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2014_I1280929);

}
