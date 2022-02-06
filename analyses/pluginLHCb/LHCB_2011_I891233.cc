// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief J/psi production at 7 TeV
  class LHCB_2011_I891233 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2011_I891233);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      vector<double> ybins={2.0,2.5,3.0,3.5,4.0,4.5};
      for(unsigned int iy=0;iy<5;++iy) {
	for(unsigned int ix=0;ix<4;++ix) {
	  Histo1DPtr tmp;
	  _h_Jpsi[ix].add(ybins[iy],ybins[iy+1],book(tmp,6+ix,1,iy+1));
	}
	Histo1DPtr tmp;
	_h_Jpsi[4].add(ybins[iy],ybins[iy+1],book(tmp,"TMP/JPsi_"+toString(iy),refData(6,1,iy+1)));
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
	if(nonPrompt) {
	  _h_Jpsi[1].fill(absrap,xp);
	}
	else {
	  _h_Jpsi[0].fill(absrap,xp);
	  _h_Jpsi[2].fill(absrap,xp);
	  _h_Jpsi[3].fill(absrap,xp);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 1/2 due rapidity folding +/-
      double factor = 0.5*crossSection()/nanobarn/sumOfWeights();
      for(unsigned int ix=0;ix<5;++ix) {
	_h_Jpsi[ix].scale(factor,this);
      }
      for(unsigned int ix=0;ix<_h_Jpsi[4].histos().size();++ix) {
	Scatter2DPtr tmp;
	book(tmp,10,1,ix);
	divide(_h_Jpsi[1].histos()[ix],_h_Jpsi[4].histos()[ix],tmp);
	tmp->scaleY(100.);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistogram _h_Jpsi[5];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2011_I891233);

}
