// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B0 production
  class CMS_2011_I896211 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2011_I896211);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_total,1,1,1);
      book(_h_pT   ,2,1,1);
      book(_h_y    ,3,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // loop over onium states
      for( const Particle & p : ufs.particles(Cuts::abspid==511)) {
	// cuts on pT and rapidity
	double y  = p.absrap();
	double pT = p.perp();
	if( pT<5. || y>2.4 ) continue;
	_h_total->fill(1.);
	_h_pT->fill(pT);
	_h_y->fill(y);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = 0.5*crossSection()/microbarn/sumOfWeights();
      scale(_h_total,fact);
      scale(_h_pT   ,fact);
      // 0.5 due +/- y
      scale(_h_y    ,0.5*fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_total,_h_pT,_h_y;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMS_2011_I896211);

}
