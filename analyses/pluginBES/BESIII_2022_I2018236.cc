// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> gamma eta' pi+pi-
  class BESIII_2022_I2018236 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2018236);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==443);
      declare(ufs, "UFS");
      DecayedParticles psi(ufs);
      psi.addStable(PID::PI0);
      psi.addStable(PID::K0S);
      psi.addStable(PID::ETA);
      psi.addStable(PID::ETAPRIME);
      declare(psi, "psi");
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  for(unsigned int iz=1;iz<2;++iz)
	    book(_h[ix][iy][iz],1+ix,1+iy,1+iz);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,1}, {-211,1}, { 331,1} ,{22,1} };
      DecayedParticles psi = apply<DecayedParticles>(event, "psi");
      // loop over particles
      for(unsigned int ix=0;ix<psi.decaying().size();++ix) {
	if(!psi.modeMatches(ix,4,mode)) continue;
	const Particle & eta  = psi.decayProducts()[ix].at( 331)[0];
	const Particle & pip  = psi.decayProducts()[ix].at( 211)[0];
	const Particle & pim  = psi.decayProducts()[ix].at(-211)[0];
	double m1 = (pip.momentum()+pim.momentum()+eta.momentum()).mass();
	double m2 = (pip.momentum()+pim.momentum()).mass();
	for(unsigned int ix=0;ix<2;++ix) {
	  _h[0][ix][1]->fill(m1);
	  _h[1][ix][1]->fill(m2);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  for(unsigned int iz=1;iz<2;++iz)
	    normalize(_h[ix][iy][iz],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2][2][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2018236);

}