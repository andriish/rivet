// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B+ meson at 13 TeV
  class CMS_2017_I1485195 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2017_I1485195);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_pT   ,5,1,3);
      book(_h_y    ,6,1,3);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // loop over B+ states
      for( const Particle & p : ufs.particles(Cuts::abspid==521)) {
	// cuts on pT and rapidity
	double y  = p.absrap();
	double pT = p.perp();
	if ( (pT<17. && y<1.45) || (pT>=17. && y<2.1) )
	  _h_pT->fill(pT);
	if ( (pT>10. && pT<100. && y<1.45) || (pT>17. && pT<100. && y>=1.45 && y<2.1)) 
	  _h_y->fill(y);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = 0.5*crossSection()/microbarn/sumOfWeights();
      scale(_h_pT   ,fact);
      // 0.5 from +/- y
      scale(_h_y    ,0.5*fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pT,_h_y;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMS_2017_I1485195);

}
