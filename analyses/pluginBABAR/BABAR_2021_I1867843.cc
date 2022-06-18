// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief eta_c Dalitz decays
  class BABAR_2021_I1867843 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2021_I1867843);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      book(_h_KK    ,1,1,1);
      book(_h_etaPK ,1,1,2);
      book(_h_pipi1 ,2,1,1);
      book(_h_etaPpi,2,1,2);
      book(_h_pipi2 ,3,1,1);
      book(_h_etapi ,3,1,2);
      book(_dalitz1, "dalitz1",50,2.,6.5,50,2.,6.5);
      book(_dalitz2, "dalitz2",50,0.,9. ,50,0.,9. );
      book(_dalitz3, "dalitz3",50,0.,8. ,50,0.,8. );
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , Particles & pi0  ,
			   Particles & Kp  , Particles & Km  , Particles & eta,
			   Particles & etaP) {
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
	else if (id == PID::ETAPRIME) {
	  etaP.push_back(p);
	  ++nstable;
	}
	else if (id == PID::K0S||id == PID::K0L) {
          ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, pi0, Kp , Km, eta,etaP);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 441 )) {
	unsigned int nstable(0);
	Particles pip, pim, pi0, Kp , Km, eta,etaP;
	findDecayProducts(meson, nstable, pip, pim, pi0, Kp , Km, eta, etaP);
	if(nstable !=3) continue;
	if(meson.pid()<0) {
	  swap(pim,pip);
	  swap(Kp,Km);
	}
	if(nstable!=3) continue;
	if (Km.size()==1&&Kp.size()==1&&etaP.size()==1&&
	     meson.mass()>2.93 && meson.mass()<3.03) {
	  double mplus  = (Kp[0].momentum()+etaP[0].momentum()).mass2();
	  double mminus = (Km[0].momentum()+etaP[0].momentum()).mass2();
	  double mKK    = (Kp[0].momentum()+Km [0].momentum()).mass2();
	  _h_KK   ->fill(sqrt(mKK));
	  _h_etaPK->fill(sqrt(mplus));
	  _h_etaPK->fill(sqrt(mminus));
	  _dalitz1->fill(mplus,mminus);
	}
	else if (pim.size()==1&&pip.size()==1&&etaP.size()==1&&
		 meson.mass()>2.93 && meson.mass()<3.03) {
	  double mplus  = (pip[0].momentum()+etaP[0].momentum()).mass2();
	  double mminus = (pim[0].momentum()+etaP[0].momentum()).mass2();
	  double mpipi    = (pip[0].momentum()+pim [0].momentum()).mass2();
	  _h_pipi1 ->fill(sqrt(mpipi));
	  _h_etaPpi->fill(sqrt(mplus));
	  _h_etaPpi->fill(sqrt(mminus));
	  _dalitz2->fill(mplus,mminus);
	}
	else if (pim.size()==1&&pip.size()==1&&eta.size()==1&&
		 meson.mass()>2.92 && meson.mass()<3.02) {
	  double mplus  = (pip[0].momentum()+eta[0].momentum()).mass2();
	  double mminus = (pim[0].momentum()+eta[0].momentum()).mass2();
	  double mpipi    = (pip[0].momentum()+pim [0].momentum()).mass2();
	  _h_pipi2->fill(sqrt(mpipi));
	  _h_etapi->fill(sqrt(mplus));
	  _h_etapi->fill(sqrt(mminus));
	  _dalitz3->fill(mplus,mminus);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_KK    );
      normalize(_h_etaPK );
      normalize(_h_pipi1 );
      normalize(_h_etaPpi);
      normalize(_h_pipi2 );
      normalize(_h_etapi );
      normalize(_dalitz1 );
      normalize(_dalitz2 );
      normalize(_dalitz3 );
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_KK,_h_etaPK;
    Histo1DPtr _h_pipi1,_h_etaPpi;
    Histo1DPtr _h_pipi2,_h_etapi;
    Histo2DPtr _dalitz1,_dalitz2,_dalitz3;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2021_I1867843);

}
