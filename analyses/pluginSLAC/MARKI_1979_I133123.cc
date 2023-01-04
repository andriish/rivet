// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D0 D+ and 7 GeV
  class MARKI_1979_I133123 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKI_1979_I133123);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::abspid==411 or Cuts::abspid==421), "UFS");
      //Histograms
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_sigma[ix],1,1,1+ix);
	for(unsigned int iy=0;iy<2;++iy)
	  book(_h_spect[ix][iy],3+iy,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
	const double xE = 2.*p.E()/sqrtS();
	unsigned int iloc = p.abspid()==421 ? 0 : 1;
	_h_sigma[iloc]->fill(sqrtS());
	for(unsigned int iy=0;iy<2;++iy)
	  _h_spect[iloc][iy]->fill(xE);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	scale(_h_sigma[ix],1.8*crossSection()/nanobarn/sumOfWeights());
	for(unsigned int iy=0;iy<2;++iy)
	  scale(_h_spect[ix][iy],sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_sigma[2],_h_spect[2][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(MARKI_1979_I133123);

}
