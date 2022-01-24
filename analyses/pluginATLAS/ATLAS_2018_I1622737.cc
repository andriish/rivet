// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief onium production at 5.02 TeV
  class ATLAS_2018_I1622737 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2018_I1622737);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      // binning in y
      vector<double> ybins={0.,0.75,1.5,2.0};
      // book histos
      for(unsigned int ix=0;ix<7;++ix) {
	for(unsigned int iy=0;iy<3;++iy) {
	  Histo1DPtr tmp;
	  _h_Onium[ix].add(ybins[iy],ybins[iy+1],book(tmp,ix+1,1,iy+1));
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles(Cuts::pid==443 or Cuts::pid==100443 or
					     Cuts::pid==553 or Cuts::pid==100553 or Cuts::pid==200553)) {
	double absrap = p.absrap();
	double xp = p.perp();
	if(p.pid()==443 || p.pid()==100443) {
	  // prompt/non-prompt
	  bool prompt = !p.fromBottom();
	  if(p.pid()==443) {
	    _h_Onium[2*prompt].fill(absrap,xp);
	  }
	  else {
	    _h_Onium[1+2*prompt].fill(absrap,xp);
	  }
	}
	else if(p.pid()==553)
	  _h_Onium[4].fill(absrap,xp);
	else if(p.pid()==100553)
	  _h_Onium[5].fill(absrap,xp);
	else if(p.pid()==200553)
	  _h_Onium[6].fill(absrap,xp);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 2 due rapidity folding +/-
      double factor = 2.*crossSection()/nanobarn/sumOfWeights();
      // br to muons PDG 2021 (psi2s is e+e- due large errors on mu+mu-)
      vector<double> br = {0.05961,0.00793,0.05961,0.00793,0.0248,0.0193,0.0218};
      for(unsigned int ix=0;ix<7;++ix) {
	_h_Onium[ix].scale(factor*br[ix],this);
      }										  
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistogram _h_Onium[7];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ATLAS_2018_I1622737);

}
