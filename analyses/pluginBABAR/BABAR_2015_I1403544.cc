// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief eta_c -> K+K-pi0 KS0 K+-pi-+
  class BABAR_2015_I1403544 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2015_I1403544);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_Kppip,3,1,1);
      book(_h_K0pip,3,1,2);
      book(_h_K0Kp ,3,1,3);
      book(_dalitz[0], "dalitz_1",50,0.,8.,50,0.,8.);
      book(_h_Kppi0,4,1,1);
      book(_h_Kmpi0,4,1,2);
      book(_h_KpKm ,4,1,3);
      book(_dalitz[1], "dalitz_2",50,0.,8.,50,0.,8.);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , Particles & pi0  ,
			   Particles & Kp  , Particles & Km  , Particles & KS0) {
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
	else if (id == PID::K0S) {
	  KS0.push_back(p);
	  ++nstable;
	}
	else if (id == PID::ETA||id == PID::K0L) {
          ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, pi0, Kp , Km, KS0);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 441 )) {
	unsigned int nstable(0);
	Particles pip, pim, pi0, Kp , Km, KS0;
	findDecayProducts(meson, nstable, pip, pim, pi0, Kp , Km, KS0);
	if(nstable !=3) continue;
	if (Km.size()==1&&Kp.size()==1&&pi0.size()==1&&
	    meson.mass()>2.922 && meson.mass()<3.036) {
	  double mplus  = (Kp[0].momentum()+pi0[0].momentum()).mass2();
	  double mminus = (Km[0].momentum()+pi0[0].momentum()).mass2();
	  double mKK    = (Kp[0].momentum()+Km [0].momentum()).mass2();
	  _h_KpKm->fill(mKK);
	  _h_Kppi0->fill(mplus);
	  _h_Kmpi0->fill(mminus);
	  _dalitz[1]->fill(mplus,mminus);
	}
	else if (((Km.size()==1&&pip.size()==1) ||
		  (Kp.size()==1&&pim.size()==1) )&&KS0.size()==1&&
		 meson.mass()>2.922 && meson.mass()<3.039) {
	  if(Km.size()==1) {
	    swap(Km,Kp);
	    swap(pim,pip);
	  }
	  double mplus  = (Kp[0].momentum()  + pim[0].momentum()).mass2();
	  double mminus = (KS0[0].momentum() + pim[0].momentum()).mass2();
	  double mKK    = (Kp[0].momentum()  + KS0[0].momentum()).mass2();
	  _h_K0Kp ->fill(mKK);
	  _h_Kppip->fill(mplus);
	  _h_K0pip->fill(mminus);
	  _dalitz[0]->fill(mplus,mminus);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Kppip);
      normalize(_h_K0pip);
      normalize(_h_K0Kp );
      normalize(_dalitz[0]);
      normalize(_h_Kppi0);
      normalize(_h_Kmpi0);
      normalize(_h_KpKm );
      normalize(_dalitz[1]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kppip,_h_K0pip,_h_K0Kp;
    Histo1DPtr _h_Kppi0,_h_Kmpi0,_h_KpKm;
    Histo2DPtr _dalitz[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2015_I1403544);

}
