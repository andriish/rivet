// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief chi_b(1p) production at 7 TeV
  class LHCB_2012_I1184177 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2012_I1184177);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_pT[ix],"TMP/pT_"+toString(ix),refData(1,1,1));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles(Cuts::pid==553 ||
					     Cuts::pid==10551 ||
					     Cuts::pid==20553 ||
					     Cuts::pid==555)) {
	// J/psi as a reference
	if(p.pid()==553) {
	  double absrap=p.absrap();
	  if(absrap>2. && absrap<4.5) _h_pT[1]->fill(p.perp());
	}
	else {
	  Particle Upsilon;
	  if(p.children()[0].pid()==22 && p.children()[1].pid()==553) {
	    Upsilon=p.children()[1];
	  }
	  else if(p.children()[1].pid()==22 && p.children()[0].pid()==553) {
	    Upsilon=p.children()[0];
	  }
	  else
	    continue;
	  double absrap=Upsilon.absrap();
	  if(absrap<2. || absrap>4.5) continue;
	  _h_pT[0]->fill(Upsilon.perp());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // compute fraction
      Scatter2DPtr tmp;
      book(tmp,1,1,1);
      efficiency(_h_pT[0],_h_pT[1],tmp);
      tmp->scaleY(100.);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pT[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2012_I1184177);

}
