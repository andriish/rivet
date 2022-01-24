// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief J/psi at 5.02 TeV
  class CMS_2017_I1512296 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2017_I1512296);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      // histograms
      vector<double> ybins={-2.4,-1.93,-1.5,-0.9,0.,0.9,1.5,1.93,2.4};
      for(int ix=0;ix<2;++ix) {
	for(int iy=0;iy<8;++iy) {
	  Histo1DPtr tmp;
	  if(iy<4)
	    _h_JPsi_pT[ix] .add(ybins[iy],ybins[iy+1],book(tmp,2+13*ix,1,4-iy));
	  else
	    _h_JPsi_pT[ix] .add(ybins[iy],ybins[iy+1],book(tmp,1+13*ix,1,8-iy));
	}
      }
      vector<double> pTbins={6.5,10.,30.};
      for(int ix=0;ix<2;++ix) {
	for(int iy=0;iy<2;++iy) {
	  Histo1DPtr tmp;
	  _h_JPsi_y[ix] .add(pTbins[iy],pTbins[iy+1],book(tmp,5+13*ix,1,iy+1));
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles(Cuts::pid==443)) {
	// prompt/non-prompt
	bool nonPrompt = p.fromBottom();
        double rap = p.rap();
        double xp = p.perp();
	_h_JPsi_pT[nonPrompt].fill(rap,xp);
	_h_JPsi_y [nonPrompt].fill(xp,rap);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // br to muons PDG 2021
      double br = 0.05961;
      double factor = br*crossSection()/microbarn/sumOfWeights();
      vector<double> pTbins={6.5,10.,30.};
      for(unsigned int ix=0;ix<2;++ix) {
	_h_JPsi_pT[ix].scale(factor,this);
	_h_JPsi_y [ix].scale(factor,this);
	// in y just single differential, multiply by bin width
	for(unsigned int iy=0;iy<2;++iy) {
	  scale(_h_JPsi_y [ix].histos()[iy],pTbins[iy+1]-pTbins[iy]);
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistogram _h_JPsi_pT[2],_h_JPsi_y[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMS_2017_I1512296);

}
