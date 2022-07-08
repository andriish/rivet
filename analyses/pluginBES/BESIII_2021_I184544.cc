// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  Ds -> KS0 K-pi+pi+
  class BESIII_2021_I184544 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I184544);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<9;++ix)
	book(_h[ix],1,1,1+ix);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , Particles & K0  ,
			   Particles & Kp  , Particles & Km  ) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::KPLUS ) {
       	  Kp.push_back(p);
	  ++nstable;
	}
	else if (id == PID::KMINUS ) {
	  Km.push_back(p);
	  ++nstable;
	}
	else if (id == PID::PIPLUS) {
	  pip.push_back(p);
	  ++nstable;
	}
	else if (id == PID::PIMINUS) {
	  pim.push_back(p);
	  ++nstable;
	}
	else if (id == PID::K0S) {
	  K0.push_back(p);
	  ++nstable;
	}
	else if (id == PID::PI0||id == PID::K0L) {
          ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, K0, Kp , Km);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 431)) {
	unsigned int nstable(0);
	Particles pip, pim, K0, Kp , Km;
	findDecayProducts(meson, nstable, pip, pim, K0, Kp , Km);
	if (nstable!=4) continue;
	if(meson.pid()<0) {
	  swap(pim,pip);
	  swap(Kp,Km);
	}
	if(pip.size()==2&&Km.size()==1&&K0.size()==1) {
	  double mK0pi[2] = {(K0[0].momentum()+pip[0].momentum()).mass(), 
			     (K0[0].momentum()+pip[1].momentum()).mass()};
	  if(mK0pi[0]>mK0pi[1]) {
	    swap(  pip[0],  pip[1]);
	    swap(mK0pi[0],mK0pi[1]);
	  }
	  _h[0]->fill((K0 [0].momentum()+Km [0].momentum()).mass());
	  _h[1]->fill(mK0pi[0]);
	  _h[2]->fill((Km [0].momentum()+pip[1].momentum()).mass());
	  _h[3]->fill(mK0pi[1]);
	  _h[4]->fill((Km [0].momentum()+pip[0].momentum()).mass());
	  _h[5]->fill((K0 [0].momentum()+Km [0].momentum()+pip[0].momentum()).mass());
	  _h[6]->fill((K0 [0].momentum()+Km [0].momentum()+pip[1].momentum()).mass());
	  _h[7]->fill((pip[0].momentum()+pip[1].momentum()).mass());
	  _h[8]->fill((Km [0].momentum()+pip[0].momentum()+pip[1].momentum()).mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<9;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[9];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I184544);

}
