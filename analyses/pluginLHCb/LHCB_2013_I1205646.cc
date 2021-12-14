// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/psi production at 2.76 TeV
  class LHCB_2013_I1205646 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2013_I1205646);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(), "UFS");
      // histos
      book(_h_pT,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // loop over onium states
      for( const Particle & p : ufs.particles(Cuts::pid==443)) {
	// cuts on pT and rapidity
	double y = p.absrap();
	double pT = p.perp();
	if(y<2. || y>4.5 || pT>12.) continue;
	_h_pT->fill(pT);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 0.5 due folded rapidity
      double factor = 0.5*crossSection() / nanobarn/ sumOfWeights();
      scale(_h_pT,factor);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_pT;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2013_I1205646);

}
