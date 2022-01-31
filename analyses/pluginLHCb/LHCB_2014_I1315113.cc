// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief chi_b production at 7 TeV
  class LHCB_2014_I1315113 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2014_I1315113);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      for(unsigned int ichi=0;ichi<2;++ichi) {
	book(_h_chi[ichi],"TMP/h_CHI_"+toString(ichi),refData(1,1,1));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles(Cuts::pid==20553 ||
					     Cuts::pid==555)) {
	// Upsilon /gamma mode
	if(p.children().size()!=2) continue;
	Particle ups;
	if(p.children()[0].pid()==22 && p.children()[1].pid()==553) {
	  ups  = p.children()[1];
	}
	else if(p.children()[1].pid()==22 && p.children()[0].pid()==553) {
	  ups  = p.children()[0];
	}
	else
	  continue;
	double absrap=ups.absrap();
	if(absrap<2. || absrap>4.5) continue;
	unsigned int ichi = p.pid()==20553 ? 0 : 1;
	_h_chi[ichi]->fill(ups.perp());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      vector<double> br = {.352,0.18};
      for(unsigned int ix=0;ix<2;++ix)
	scale(_h_chi[ix],1./br[ix]);
      Scatter2DPtr tmp;
      book(tmp,1,1,1);
      divide(_h_chi[0],_h_chi[1],tmp);
      tmp->scaleY(br[1]/br[0]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_chi[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2014_I1315113);

}
