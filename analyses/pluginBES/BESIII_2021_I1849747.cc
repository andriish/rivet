// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Ds -> K+K-pi+pi0
  class BESIII_2021_I1849747 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1849747);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<10;++ix)
	book(_h[ix],1,1,1+ix);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , Particles & pi0  ,
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
	else if (id == PID::PI0) {
	  pi0.push_back(p);
	  ++nstable;
	}
	else if (id == PID::K0S||id == PID::K0L) {
          ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, pi0, Kp , Km);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 431)) {
	unsigned int nstable(0);
	Particles pip, pim, pi0, Kp , Km;
	findDecayProducts(meson, nstable, pip, pim, pi0, Kp , Km);
	if (nstable!=4) continue;
	if(meson.pid()<0) {
	  swap(pim,pip);
	  swap(Kp,Km);
	}
	if(pip.size()==1&&Km.size()==1&&Kp.size()==1&&pi0.size()==1) {
	  double mKK = (Kp[0].momentum()+Km[0].momentum()).mass();
	  _h[0]->fill(mKK);
	  _h[1]->fill(mKK);
	  _h[2]->fill((Kp [0].momentum()+pi0[0].momentum()).mass());
	  _h[3]->fill((Km [0].momentum()+pi0[0].momentum()).mass());
	  _h[4]->fill((pip[0].momentum()+pi0[0].momentum()).mass());
	  _h[5]->fill((Km [0].momentum()+pip[0].momentum()).mass());
	  _h[6]->fill((Km [0].momentum()+pip[0].momentum()+pi0[0].momentum()).mass());
	  _h[7]->fill((Km [0].momentum()+ Kp[0].momentum()+pip[0].momentum()).mass());
	  _h[8]->fill((Kp [0].momentum()+pip[0].momentum()+pi0[0].momentum()).mass());
	  _h[9]->fill((Km [0].momentum()+ Kp[0].momentum()+pi0[0].momentum()).mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<10;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[10];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1849747);

}
