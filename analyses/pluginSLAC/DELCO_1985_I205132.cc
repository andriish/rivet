// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D* +/- production at 29 GeV
  class DELCO_1985_I205132 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DELCO_1985_I205132);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::abspid==413), "UFS");
      //Histograms
      for(unsigned int ix=0;ix<2;++ix) book(_h[ix],1,1,1+ix);
      book(_r,2,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
	const double xp = p.p3().mod()/sqrt(0.25*sqr(sqrtS())-p.mass2());
	if(xp>0.3) _r->fill(sqrtS());
	_h[0]->fill(xp);
	if(!p.fromBottom()) _h[1]->fill(xp);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	scale(_h[ix], sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
      scale(_r,crossSection()/nanobarn/sumOfWeights());
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2], _r;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(DELCO_1985_I205132);

}
