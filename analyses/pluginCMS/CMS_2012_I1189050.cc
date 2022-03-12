// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief chi_c at 7 TeV
  class CMS_2012_I1189050 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2012_I1189050);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      for(unsigned int ichi=0;ichi<2;++ichi) {
	for(unsigned int ipT=0;ipT<2;++ipT) {
	  book(_h_chi[ichi][ipT],"TMP/h_CHI_"+toString(ichi)+"_"+toString(ipT),refData(1,1,1));
	}
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
	Particle Jpsi,gamma;
	if(p.children()[0].pid()==22 && p.children()[1].pid()==443) {
	  Jpsi  = p.children()[1];
	  gamma = p.children()[0];
	}
	else if(p.children()[1].pid()==22 && p.children()[0].pid()==443) {
	  Jpsi  = p.children()[0];
	  gamma = p.children()[1];
	}
	else
	  continue;
	if(Jpsi.absrap()>1.) continue;
	double xp=Jpsi.perp();
	unsigned int ichi = p.pid()==20443 ? 0 : 1;
	if(gamma.perp()>0.5)
	  _h_chi[ichi][0]->fill(xp);
	_h_chi[ichi][1]->fill(xp);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // chi_c to gamma J/psi branching ratios from PDG 2021
      vector<double> br = {0.343,0.190};
      // no br correction fiduical region
      Scatter2DPtr tmp;
      book(tmp,1,1,1);
      divide(_h_chi[1][0],_h_chi[0][0],tmp);
      // br correction fiduical region
      book(tmp,2,1,1);
      divide(_h_chi[1][0],_h_chi[0][0],tmp);
      tmp->scaleY(br[0]/br[1]);
      // no br correction  pT gamma >0
      book(tmp,3,1,1);
      divide(_h_chi[1][1],_h_chi[0][1],tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_chi[2][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMS_2012_I1189050);

}
