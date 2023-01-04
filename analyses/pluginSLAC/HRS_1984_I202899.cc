// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D0 and D+ at 29 GeV
  class HRS_1984_I202899 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(HRS_1984_I202899);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::abspid==411 or Cuts::abspid==421), "UFS");
      //Histograms
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_sigma[ix],1,1,1+ix);
	for(unsigned int iy=0;iy<2;++iy) {
	  if(ix==1&&iy==1) continue;
	  book(_h_spect[ix][iy],3+iy,1,1+ix);
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
       	const double xE = 2.*p.E()/sqrtS();
       	const double beta = p.p3().mod() / p.E();
      	unsigned int iloc = p.abspid()==421 ? 0 : 1;
       	if(xE>0.5) _h_sigma[iloc]->fill(sqrtS());
	_h_spect[iloc][0]->fill(xE,1./beta);
	if(iloc==0) _h_spect[0][1]->fill(xE);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	scale(_h_sigma[ix],crossSection()/nanobarn/sumOfWeights());
	scale(_h_spect[ix][0],sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
      }
      normalize(_h_spect[0][1]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_sigma[2],_h_spect[2][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(HRS_1984_I202899);

}
