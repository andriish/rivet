// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief  J/psi production at 13 TeV
  class LHCB_2015_I1391511 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2015_I1391511);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      vector<double> ybins={2.0,2.5,3.0,3.5,4.0,4.5};
      for(unsigned int iy=0;iy<5;++iy) {
	for(unsigned int ix=0;ix<2;++ix) {
	  Histo1DPtr tmp;
	  _h_Jpsi[ix].add(ybins[iy],ybins[iy+1],book(tmp,1+ix,1,iy+1));
	}
	Histo1DPtr tmp;
	_h_Jpsi[2].add(ybins[iy],ybins[iy+1],book(tmp,"TMP/Jpsi_"+toString(iy),refData(3,1,iy+1)));
      }
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_pT[ix],4+ix,1,1);
	book(_h_y [ix],6+ix,1,1);
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
	_h_Jpsi[nonPrompt].fill(absrap,xp);
	_h_Jpsi[2      ].fill(absrap,xp);
	_h_pT[nonPrompt]->fill(xp);
	_h_y [nonPrompt]->fill(absrap);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 1/2 due rapidity folding +/-
      double factor = 0.5*crossSection()/nanobarn/sumOfWeights();
      for(unsigned int ix=0;ix<3;++ix) {
	_h_Jpsi[ix].scale(factor,this);
      }
      Scatter2DPtr tmp;
      for(unsigned int ix=0;ix<_h_Jpsi[1].histos().size();++ix) {
	book(tmp,3,1,1+ix);
	divide(_h_Jpsi[1].histos()[ix],_h_Jpsi[2].histos()[ix],tmp);
	tmp->scaleY(100.);
      }
      for(unsigned int ix=0;ix<2;++ix) {
	scale(_h_pT[ix],factor);
	scale(_h_y [ix],factor*1e-3);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistogram _h_Jpsi[3];
    Histo1DPtr _h_pT[2],_h_y[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2015_I1391511);

}
