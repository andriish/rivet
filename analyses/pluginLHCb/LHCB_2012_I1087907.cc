// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief chi_c at 7 TeV
  class LHCB_2012_I1087907 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2012_I1087907);


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

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles(Cuts::pid==20443 ||
					     Cuts::pid==445)) {
	// prompt
	if(p.fromBottom()) continue;
	// J/psi /gamma mode
	if(p.children().size()!=2) continue;
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
	unsigned int ichi = p.pid()==20443 ? 0 : 1;
	double xp=Jpsi.perp();
	_h_chi[ichi]->fill(xp);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // chi_c to gamma J/psi branching ratios from PDG 2021
      vector<double> br = {0.343,0.190};
      Scatter2DPtr tmp;
      book(tmp,1,1,1);
      divide(_h_chi[1],_h_chi[0],tmp);
      tmp->scaleY(br[0]/br[1]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_chi[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2012_I1087907);

}
