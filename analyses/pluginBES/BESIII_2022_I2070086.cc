// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D_s+ -> KS0 K+ pi0
  class BESIII_2022_I2070086 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2070086);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_KK  ,1,1,1);
      book(_h_K0pi,1,1,2);
      book(_h_Kppi,1,1,3);
      book(_dalitz,"dalitz",50,0.3,2.3,50,0.3,2.3);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & Kp , Particles & Km , Particles & pi0,
			   Particles & K0) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::PIPLUS || id == PID::PIMINUS ||
	     id == PID::ETA   || id == PID::K0L) {
	  ++nstable;
	}
	else if (id == PID::KPLUS) {
	  Kp.push_back(p);
	  ++nstable;
	}
	else if (id == PID::KMINUS) {
	  Km.push_back(p);
	  ++nstable;
	}
	else if (id == PID::PI0) {
	  pi0.push_back(p);
	  ++nstable;
	}
	else if (id == PID::K0S) {
	  K0.push_back(p);
          ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, Kp, Km, pi0, K0);
	}
	else
	  ++nstable;
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 431 )) {
      	unsigned int nstable(0);
      	Particles Kp, Km, pi0, K0;
      	findDecayProducts(meson, nstable, Kp, Km, pi0, K0);
      	if(nstable !=3) continue;
      	if(meson.pid()<0) {
      	  swap(Km,Kp);
      	}
      	if (K0.size()==1&&Kp.size()==1&&pi0.size()==1) {
     	  double m0  = (K0[0].momentum()+pi0[0].momentum()).mass2();
     	  double mp  = (Kp[0].momentum()+pi0[0].momentum()).mass2();
     	  double mKK = (K0[0].momentum()+Kp[0].momentum()).mass2();
      	  _dalitz->fill(mp,m0);
      	  _h_KK->fill(sqrt(mKK));
	  _h_Kppi->fill(sqrt(mp));
	  _h_K0pi->fill(sqrt(m0));
      	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_K0pi);
      normalize(_h_Kppi);
      normalize(_h_KK  );
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_K0pi,_h_Kppi,_h_KK;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2070086);

}
