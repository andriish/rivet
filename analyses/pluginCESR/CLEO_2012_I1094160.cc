// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D0 -> KS) K+/- pi-/+
  class CLEO_2012_I1094160 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2012_I1094160);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // histograms
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_K0Km [ix],1,1,1+3*ix);
	book(_h_K0pip[ix],1,1,2+3*ix);
	book(_h_Kmpip[ix],1,1,3+3*ix);
	book(_h_K0Kp [ix],2,1,1+3*ix);
	book(_h_K0pim[ix],2,1,2+3*ix);
	book(_h_Kppim[ix],2,1,3+3*ix);
	book(_dalitz [ix],"dalitz_"+toString(ix+1),50,0.3,2.0,50,0.3,2.);
      }
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , Particles & K0, Particles & Kp,  Particles & Km) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if (id == PID::PI0 || id == PID::K0L) {
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
	else if (id == PID::KPLUS) {
	  Kp.push_back(p);
	  ++nstable;
	}
	else if (id == PID::KMINUS) {
	  Km.push_back(p);
	  ++nstable;
	}
	else if (id == PID::K0S) {
	  K0.push_back(p);
          ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, K0, Kp, Km);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").
	    particles(Cuts::abspid== 421)) {
	unsigned int nstable(0);
	Particles pip, pim, K0,Kp,Km;
	findDecayProducts(meson, nstable, pip, pim, K0,Kp,Km);
	if(nstable !=3 || K0.size()!=1) continue;
	if(meson.pid()<0) {
	  swap(pim,pip);
	  swap(Km,Kp);
	}
	// K0S K- pi+
	if( Km.size()==1 && pip.size()==1) {
	  double mK0pip = (K0[0].momentum()+pip[0].momentum() ).mass2();
	  double mKmpip = (Km[0].momentum()+pip[0].momentum() ).mass2();
	  double mKK    = (K0[0].momentum()+Km [0].momentum() ).mass2();
	  for(unsigned int ix=0;ix<2;++ix) {
	    _h_K0Km [ix]->fill(mKK   );
	    _h_K0pip[ix]->fill(mK0pip);
	    _h_Kmpip[ix]->fill(mKmpip);
	  }
	  _dalitz[0]->fill(mKmpip,mK0pip); 
	}
	// K0S K+ pi-
	else if( Kp.size()==1 && pim.size()==1) {
	  double mK0pim = (K0[0].momentum()+pim[0].momentum() ).mass2();
	  double mKppim = (Kp[0].momentum()+pim[0].momentum() ).mass2();
	  double mKK    = (K0[0].momentum()+Kp [0].momentum() ).mass2();
	  for(unsigned int ix=0;ix<2;++ix) {
	    _h_K0Kp [ix]->fill(mKK   );
	    _h_K0pim[ix]->fill(mK0pim);
	    _h_Kppim[ix]->fill(mKppim);
	  }
	  _dalitz[0]->fill(mKppim,mK0pim); 
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_h_K0Km [ix]);
	normalize(_h_K0pip[ix]);
	normalize(_h_Kmpip[ix]);
	normalize(_h_K0Kp [ix]);
	normalize(_h_K0pim[ix]);
	normalize(_h_Kppim[ix]);
	normalize(_dalitz [ix]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kmpip[2], _h_K0pip[2], _h_K0Km[2];
    Histo1DPtr _h_Kppim[2], _h_K0pim[2], _h_K0Kp[2];
    Histo2DPtr _dalitz[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_2012_I1094160);

}
