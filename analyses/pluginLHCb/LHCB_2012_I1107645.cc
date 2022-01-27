// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief chi_c production at 7 TeV
  class LHCB_2012_I1107645 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2012_I1107645);


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

      for (const Particle& p : ufs.particles(Cuts::pid==443 ||
					     Cuts::pid==10441 ||
					     Cuts::pid==20443 ||
					     Cuts::pid==445)) {
	// prompt only
	if(p.fromBottom()) continue;
	// J/psi as a reference
	if(p.pid()==443) {
	  double absrap=p.absrap();
	  if(absrap>2. && absrap<4.5) _h_pT[1]->fill(p.perp());
	}
	else {
	  Particle Jpsi;
	  if(p.children()[0].pid()==22 && p.children()[1].pid()==443) {
	    Jpsi=p.children()[1];
	  }
	  else if(p.children()[1].pid()==22 && p.children()[0].pid()==443) {
	    Jpsi=p.children()[0];
	  }
	  else
	    continue;
	  double absrap=Jpsi.absrap();
	  if(absrap<2. || absrap>4.5) continue;
	  _h_pT[0]->fill(Jpsi.perp());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // compute fraction
      Scatter2DPtr tmp;
      book(tmp,1,1,1);
      efficiency(_h_pT[0],_h_pT[1],tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pT[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2012_I1107645);

}
