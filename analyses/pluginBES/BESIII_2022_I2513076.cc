// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief pi0/K0 spectra
  class BESIII_2022_I2513076 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2513076);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(Cuts::pid==PID::PI0 or
				Cuts::pid==PID::K0S), "UFS");
      // find beam energy
      unsigned int iloc(0);
      if(isCompatibleWithSqrtS(2.2324, 1E-3))
	iloc = 1;
      else if(isCompatibleWithSqrtS(2.4 , 1E-3))
	iloc = 2;
      else if(isCompatibleWithSqrtS(2.8 , 1E-3))
	iloc = 3;
      else if(isCompatibleWithSqrtS(3.05 , 1E-3))
	iloc = 4;
      else if(isCompatibleWithSqrtS(3.4 , 1E-3))
	iloc = 5;
      else if(isCompatibleWithSqrtS(3.671 , 1E-3))
	iloc = 6;
      else
	MSG_ERROR("Beam energy not supported!");
      // book histos
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1+ix,1,iloc);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
	double pp = p.p3().mod();
	if(p.pid()==PID::PI0) _h[0]->fill(pp);
	else                  _h[1]->fill(pp);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	scale(_h[ix],1./sumOfWeights());
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2513076);

}
