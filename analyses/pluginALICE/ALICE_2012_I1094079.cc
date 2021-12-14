// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/psi production at 2.76 TeV
  class ALICE_2012_I1094079 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2012_I1094079);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(), "UFS");
      // histos
      book(_h_pT,1,1,1);
      book(_h_y ,4,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // loop over onium states
      for( const Particle & p : ufs.particles(Cuts::pid==443)) {
	// cuts on pT and rapidity
	double y = p.absrap();
	if(y>2.5 && y<4.)
	  _h_pT->fill(p.perp());
	_h_y->fill(p.rapidity());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double factor = crossSection() / microbarn/ sumOfWeights();
      scale(_h_pT,factor/3.);
      scale(_h_y ,factor);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_pT,_h_y;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(ALICE_2012_I1094079);

}
