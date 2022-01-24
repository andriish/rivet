// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief 
  class LHCB_2021_I1915030 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2021_I1915030);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      // histograms
      vector<double> ybins={2.0,2.5,3.0,3.5,4.0,4.5};
      for(int iy=0;iy<5;++iy) {
	for(int ix=0;ix<2;++ix) {
	  Histo1DPtr tmp;
	  _h_JPsi[ix] .add(ybins[iy],ybins[iy+1],book(tmp,ix+1,1,iy+1));
	}
	Histo1DPtr tmp;
	_h_JPsi[2].add(ybins[iy],ybins[iy+1],book(tmp ,"TMP/JPsi_"+toString(iy) , refData(1,1,iy+1)));
      }
      book(_h_JPsi_pT[0],3,1,1);
      book(_h_JPsi_pT[1],4,1,1);
      book(_h_JPsi_y [0],5,1,1);
      book(_h_JPsi_y [2],5,1,2);
      book(_h_JPsi_y [1],6,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles(Cuts::pid==443)) {
	// prompt/non-prompt
	bool nonPrompt = p.fromBottom();
        double absrap = p.absrap();
        double xp = p.perp();
	_h_JPsi[nonPrompt].fill(absrap,xp);
	_h_JPsi[2        ].fill(absrap,xp);
	if(absrap>2. && absrap<4.5)
	  _h_JPsi_pT[nonPrompt]->fill(xp);
	if(xp<14.)
	  _h_JPsi_y [nonPrompt]->fill(absrap);
	if(xp<8. && !nonPrompt)
	  _h_JPsi_y [2]->fill(absrap);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // br to muons PDG 2021
      double br = 0.05961;
      // factor of 2 due rapidity +/-
      double factor = 2.*br*crossSection()/nanobarn/sumOfWeights();
      for(unsigned int ix=0;ix<3;++ix) {
	_h_JPsi[ix].scale(factor,this);
	scale(_h_JPsi_y[ix],factor);
	if(ix<2)
	  scale(_h_JPsi_pT[ix],factor);
      }
      for(unsigned int iy=0;iy<_h_JPsi[0].histos().size();++iy) {
	// non-prompt J/psi percentage
	Scatter2DPtr tmp;
	book(tmp,7,1,iy+1);
	efficiency(_h_JPsi[1].histos()[iy],_h_JPsi[2].histos()[iy],tmp);
	tmp->scaleY(100.);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistogram _h_JPsi[3];
    Histo1DPtr _h_JPsi_pT[2],_h_JPsi_y[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2021_I1915030);

}
