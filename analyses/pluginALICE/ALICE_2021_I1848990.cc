// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D meson production at 5.02 TeV
  class ALICE_2021_I1848990 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2021_I1848990);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(), "UFS");
      // histograms
      for(unsigned int ix=0;ix<3;++ix) {
	book(_h_prompt[ix][1],1+ix,1,1);
	book(_h_prompt[ix][0],4+ix,1,1);
	book(_h_prompt[ix][2], "TMP/h_prompt"+toString(ix+1),refData(7+ix,1,1));
      }
      book(_h_D0  [0],"TMP/h_D0_2"  ,refData(10,1,1));
      book(_h_D0  [1],"TMP/h_D0_1"  ,refData(11,1,1));
      book(_h_Dsum[0],"TMP/h_Dsum_1",refData(12,1,1));
      book(_h_Dsum[1],"TMP/h_Dsum_2",refData(13,1,1));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==411 or Cuts::abspid==421 or
					     Cuts::abspid==431)) {
	// no mixing and |y|<0.5
	if(p.children().size()==1 || p.absrap()>0.5) continue;
	unsigned int iprompt = p.fromBottom();
	unsigned int imeson=0;
	if     (p.abspid()==411) imeson=1;
	else if(p.abspid()==431) imeson=2;
	double pT=p.perp();
	_h_prompt[imeson][iprompt]->fill(pT);
	if(iprompt==0)
	  _h_prompt[imeson][2]->fill(pT);
	if(imeson==0) _h_D0  [iprompt]->fill(pT);
	if(imeson <2) _h_Dsum[iprompt]->fill(pT);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double factor = crossSection()/microbarn/sumOfWeights();
      for(unsigned int ix=0;ix<3;++ix) {
      	for(unsigned int iy=0;iy<3;++iy)
      	  scale(_h_prompt[ix][iy],factor);
      }
      for(unsigned int ix=0;ix<2;++ix) {
      	scale(_h_D0  [ix],factor);
      	scale(_h_Dsum[ix],factor);
      }
      // non-prompt/prompt ratios
      for(unsigned int ix=0;ix<3;++ix) {
      	Scatter2DPtr tmp;
      	book(tmp,7+ix,1,1);
      	divide(_h_prompt[ix][1],_h_prompt[ix][2],tmp);
      }
      for(unsigned int ix=0;ix<2;++ix) {
      	Scatter2DPtr tmp;
      	// prompt and non-prompt D+/D0
      	book(tmp,10+ix,1,1);
      	divide(_h_prompt[1][ix],_h_D0[ix],tmp);
      	// prompt and non-prompt D_s/(D0+D+)
      	book(tmp,12+ix,1,1);
      	divide(_h_prompt[2][ix],_h_Dsum[ix],tmp);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_prompt[3][3];
    Histo1DPtr _h_D0[2],_h_Dsum[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ALICE_2021_I1848990);

}
