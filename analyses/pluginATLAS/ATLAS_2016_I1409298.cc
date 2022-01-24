// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief ATLAS J/psi psi2s at 7 and 8 TeV
  class ATLAS_2016_I1409298 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2016_I1409298);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      unsigned int iloc=0;
      if (isCompatibleWithSqrtS(7000)) {
	iloc = 1;
      }
      else if  (isCompatibleWithSqrtS(8000)) {
	iloc = 2;
      }
      else
	throw UserError("Centre-of-mass energy of the given input is neither 7 or 8 TeV.");
      // binning in y
      vector<double> ybins={0.,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.};
      // book histos
      for(unsigned int iy=0;iy<8;++iy) {
	for(unsigned int ix=0;ix<2;++ix) {
	  Histo1DPtr tmp;
	  // prompt and non-prompt Jpsi
	  _h_JPsi[ix] .add(ybins[iy],ybins[iy+1],book(tmp,  iloc+2*ix,1,iy+1));
	  // prompt and non-prompt psi(2S)
	  _h_psi2S[ix].add(ybins[iy],ybins[iy+1],book(tmp,4+iloc+2*ix,1,iy+1));
	}
	// total no for ratios etc
	Histo1DPtr tmp;
	_h_JPsi[2] .add(ybins[iy],ybins[iy+1],book(tmp ,"TMP/JPsi_"+toString(iy) , refData(  iloc,1,iy+1)));
	Histo1DPtr tmp2;
	_h_psi2S[2].add(ybins[iy],ybins[iy+1],book(tmp2,"TMP/psi2S_"+toString(iy), refData(4+iloc,1,iy+1)));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles(Cuts::pid==443 or Cuts::pid==100443)) {
	// prompt/non-prompt
	bool nonPrompt = p.fromBottom();
        double absrap = p.absrap();
        double xp = p.perp();
	if(p.pid()==443) {
	  _h_JPsi[nonPrompt].fill(absrap,xp);
	  _h_JPsi[2        ].fill(absrap,xp);
	}
	else {
	  _h_psi2S[nonPrompt].fill(absrap,xp);
	  _h_psi2S[2        ].fill(absrap,xp);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 2 due rapidity folding +/-
      double factor = 2.*crossSection()/nanobarn/sumOfWeights();
      // br to muons PDG 2021 (psi2s is e+e- due large errors on mu+mu-)
      vector<double> br = {0.05961,0.00793};
      // scale histos
      for(unsigned int ix=0;ix<2;++ix) {
	_h_JPsi [ix].scale(factor*br[0],this);
	_h_psi2S[ix].scale(factor*br[1],this);
      }
      // ratios, first find CMS energy
      unsigned int iloc=0;
      if (isCompatibleWithSqrtS(7000)) {
	iloc = 1;
      }
      else if  (isCompatibleWithSqrtS(8000)) {
	iloc = 2;
      }
      else
	throw UserError("Centre-of-mass energy of the given input is neither 7 or 8 TeV.");
      
      for(unsigned int iy=0;iy<_h_JPsi[0].histos().size();++iy) {
	// non-prompt J/psi percentage
	Scatter2DPtr tmp;
	book(tmp,8+iloc,1,iy+1);
	efficiency(_h_JPsi[1].histos()[iy],_h_JPsi[2].histos()[iy],tmp);
	tmp->scaleY(100.);
	// non-prompt psi2S percentage
	Scatter2DPtr tmp2;
	book(tmp2,10+iloc,1,iy+1);
	efficiency(_h_psi2S[1].histos()[iy],_h_psi2S[2].histos()[iy],tmp2);
	tmp2->scaleY(100.);
	// prompt psi(2s)/J/psi percentage
	Scatter2DPtr tmp3;
	book(tmp3,12+iloc,1,iy+1);
	divide(_h_psi2S[0].histos()[iy],_h_JPsi[0].histos()[iy],tmp3);
	tmp3->scaleY(100.);
	// non-prompt psi(2s)/J/psi percentage
	Scatter2DPtr tmp4;
	book(tmp4,14+iloc,1,iy+1);
	divide(_h_psi2S[1].histos()[iy],_h_JPsi[1].histos()[iy],tmp4);
	tmp4->scaleY(100.);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistogram _h_JPsi[3],_h_psi2S[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ATLAS_2016_I1409298);

}
