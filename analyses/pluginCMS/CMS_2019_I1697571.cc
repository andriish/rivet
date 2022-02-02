// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B_s0 at 5.02 TeV
  class CMS_2019_I1697571 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2019_I1697571);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_pT[0],1,1,1);
      book(_h_pT[1],2,1,1);
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
	if(pT<7. || pT>50.) continue;
	_h_pT[0]->fill(pT);
	_h_pT[1]->fill(pT);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 0.5 from particle/antiparticle
      double fact = 0.5*crossSection()/picobarn/sumOfWeights();
      for(unsigned int ix=0;ix<2;++ix) scale(_h_pT[ix],fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pT[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMS_2019_I1697571);

}
