// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B_s production at 7 TeV
  class CMS_2011_I914325 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2011_I914325);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_sigma,1,1,1);
      book(_h_pT   ,2,1,1);
      book(_h_y    ,3,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // loop over onium states
      for( const Particle & p : ufs.particles(Cuts::abspid==531)) {
	// skip copies due mixing
	if(p.children().size()==1 && p.children()[0].abspid()==531) continue;
	// rapidity and pT
	double y = p.absrap();
	if(y>2.4) continue;
	double pT = p.perp();
	if(pT<8. || pT>50.) continue;
	_h_sigma->fill(29.);
	_h_pT->fill(pT);
	_h_y->fill(y);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // br for J/psi phi (PDG 2021)
      double br = 1.08e-3;
      // 0.5 from particle/antiparticle
      double fact = 0.5*br*crossSection()/nanobarn/sumOfWeights();
      // remove the bin width as total cross section
      scale(_h_sigma,fact*42.);
      // normal differential cross sections
      scale(_h_pT,fact);
      scale(_h_y,fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_sigma,_h_pT,_h_y;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMS_2011_I914325);

}
