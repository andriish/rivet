// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  eta_c -> K+K- eta or pi0
  class BABAR_2014_I1287632 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2014_I1287632);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      book(_h_KK[0],1,1,1);
      book(_h_KK[1],2,1,1);
      book(_h_Kpeta,1,1,2);
      book(_h_Kmeta,1,1,3);
      book(_h_Kppi ,2,1,2);
      book(_h_Kmpi ,2,1,3);
      book(_dalitz[0], "dalitz_1",50,0.5,7.0,50,0.5 ,7.0);
      book(_dalitz[1], "dalitz_2",50,0.2,6.5,50,0.2,6.5);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , Particles & pi0  ,
			   Particles & Kp  , Particles & Km  , Particles & eta) {
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
	else if (id == PID::ETA) {
	  eta.push_back(p);
	  ++nstable;
	}
	else if (id == PID::K0S||id == PID::K0L) {
          ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, pi0, Kp , Km, eta);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 441 )) {
	unsigned int nstable(0);
	Particles pip, pim, pi0, Kp , Km, eta;
	findDecayProducts(meson, nstable, pip, pim, pi0, Kp , Km, eta);
	if(nstable !=3) continue;
	if(meson.pid()<0) {
	  swap(pim,pip);
	  swap(Kp,Km);
	}
	if(nstable!=3) continue;
	if (Km.size()==1&&Kp.size()==1&&pi0.size()==1&&
	    meson.mass()>2.922 && meson.mass()<3.036) {
	  double mplus  = (Kp[0].momentum()+pi0[0].momentum()).mass2();
	  double mminus = (Km[0].momentum()+pi0[0].momentum()).mass2();
	  double mKK    = (Kp[0].momentum()+Km [0].momentum()).mass2();
	  _h_KK[1]->fill(mKK);
	  _h_Kppi->fill(mplus);
	  _h_Kmpi->fill(mminus);
	  _dalitz[1]->fill(mplus,mminus);
	}
	else if (Km.size()==1&&Kp.size()==1&&eta.size()==1&&
		 meson.mass()>2.910 && meson.mass()<3.03) {
	  double mplus  = (Kp[0].momentum()+eta[0].momentum()).mass2();
	  double mminus = (Km[0].momentum()+eta[0].momentum()).mass2();
	  double mKK    = (Kp[0].momentum()+Km [0].momentum()).mass2();
	  _h_KK[0]->fill(mKK);
	  _h_Kpeta->fill(mplus);
	  _h_Kmeta->fill(mminus);
	  _dalitz[0]->fill(mplus,mminus);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_KK[0]);
      normalize(_h_KK[1]);
      normalize(_h_Kmpi);
      normalize(_h_Kppi);
      normalize(_h_Kmeta);
      normalize(_h_Kpeta);
      normalize(_dalitz[0]);
      normalize(_dalitz[1]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_KK[2],_h_Kmeta,_h_Kpeta,_h_Kmpi,_h_Kppi;
    Histo2DPtr _dalitz[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2014_I1287632);

}
