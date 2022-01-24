// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief psi(2S) at 5.02
  class CMS_2019_I1672011 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2019_I1672011);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      // histograms
      vector<double> ybins={0.0,0.9,1.5,1.93,2.4};
      for(int iy=0;iy<4;++iy) {
	Histo1DPtr tmp;
	_h_psi2S.add(ybins[iy],ybins[iy+1],book(tmp,3,1,4-iy));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles(Cuts::pid==100443)) {
	// prompt
	if(p.fromBottom()) continue;
        double absrap = p.absrap();
        double xp = p.perp();
	_h_psi2S.fill(absrap,xp);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // br to muons PDG 2021 (e+e- due large errors on mu+mu-)
      double br = 0.00793;
      // factor of 1/2 as +/- y
      double factor = 0.5*br*crossSection()/nanobarn/sumOfWeights();
      _h_psi2S.scale(factor,this);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistogram _h_psi2S;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMS_2019_I1672011);

}
