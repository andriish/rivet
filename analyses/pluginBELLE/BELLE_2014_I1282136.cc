// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief taU -> pi- KS0 KS0 pi0 nu
  class BELLE_2014_I1282136 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2014_I1282136);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==15);
      declare(ufs, "UFS");
      DecayedParticles TAU(ufs);
      TAU.addStable(310);
      TAU.addStable(111);
      declare(TAU, "TAU");
      // histograms
      for(unsigned int ix=0;ix<7;++ix) {
	if(ix<2) book(_h[ix],1,1,1+ix);
	else     book(_h[ix],2,1,ix-1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 310,2},{ 111,1}, {-211,1},{ 16,1}};
      static const map<PdgId,unsigned int> & modeCC = { { 310,2},{ 111,1}, { 211,1},{-16,1}};
      DecayedParticles TAU = apply<DecayedParticles>(event, "TAU");
      // loop over particles
      for(unsigned int ix=0;ix<TAU.decaying().size();++ix) {
	int sign = 1;
	if(TAU.decaying()[ix].pid()>0 &&
	   TAU.modeMatches(ix,5,mode)) sign=1;
	else if (TAU.decaying()[ix].pid()<0 &&
		 TAU.modeMatches(ix,5,modeCC)) sign=-1;
	else
	  continue;
	const Particle  & pi0 = TAU.decayProducts()[ix].at(111)[0];
	const Particle  & pim = TAU.decayProducts()[ix].at(-211*sign)[0];
	const Particles & K0  = TAU.decayProducts()[ix].at( 310);
	_h[0]->fill((pi0.momentum()+K0[0].momentum()+K0[1].momentum()).mass());
	_h[2]->fill((pim.momentum()+pi0.momentum()).mass());
	_h[4]->fill((K0[0].momentum()+K0[1].momentum()).mass());
	_h[5]->fill((pim.momentum()+K0[0].momentum()+K0[1].momentum()).mass());
	for(unsigned int ix=0;ix<2;++ix) {
	  _h[1]->fill((pim.momentum()+K0[ix].momentum()).mass());
	  _h[3]->fill((pi0.momentum()+K0[ix].momentum()).mass());
	  _h[6]->fill((pi0.momentum()+pim.momentum()+K0[ix].momentum()).mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<7;++ix) {
	normalize(_h[ix],1.,false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[7];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2014_I1282136);

}
