// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief eta_c dalitz plots
  class MC_EtaC_Dalitz : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MC_EtaC_Dalitz);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // pi+pi- eta
      book(_h1_pippim,"h01_pippim",100,0.2,2.5);
      book(_h1_pipeta,"h01_pipeta",100,0.5,3.0);
      book(_h1_pimeta,"h01_pimeta",100,0.5,3.0);
      book(_dalitz[0], "dalitz01",50,0.,9.,50,0.,9.);
      // pi0 pi0 eta
      book(_h2_pi0pi0,"h02_pi0pi0",100,0.2,2.5);
      book(_h2_pi0eta,"h02_pi0eta",100,0.5,3.0);
      book(_dalitz[1], "dalitz02",50,0.,9.,50,0.,9.);
      // pi+pi- eta'
      book(_h3_pippim,"h03_pippim",100,0.2,2.2);
      book(_h3_pipeta,"h03_pipeta",100,1.0,3.0);
      book(_h3_pimeta,"h03_pimeta",100,1.0,3.0);
      book(_dalitz[2], "dalitz03",50,0.,9.,50,0.,9.);
      // pi0 pi0 eta'
      book(_h4_pi0pi0,"h04_pi0pi0",100,0.2,2.2);
      book(_h4_pi0eta,"h04_pi0eta",100,1.0,3.0);
      book(_dalitz[3], "dalitz04",50,0.,9.,50,0.,9.);
      // K+K- eta
      book(_h5_KpKm ,"h05_KpKm" ,100,0.5,2.5);
      book(_h5_Kpeta,"h05_Kpeta",100,1.0,3.0);
      book(_h5_Kmeta,"h05_Kmeta",100,1.0,3.0);
      book(_dalitz[4], "dalitz05",50,2.,7.,50,2.,7.);
      // KS0 KS0 eta
      book(_h6_KS0KS0,"h06_KS0KS0",100,0.5,2.5);
      book(_h6_KS0eta,"h06_KS0eta",100,1.0,3.0);
      book(_dalitz[5], "dalitz06",50,2.,7.,50,2.,7.);
      // KL0 KL0 eta
      book(_h7_KL0KL0,"h07_KL0KL0",100,0.5,2.5);
      book(_h7_KL0eta,"h07_KL0eta",100,1.0,3.0);
      book(_dalitz[6], "dalitz07",50,2.,7.,50,2.,7.);
      // K+K- eta'
      book(_h8_KpKm ,"h08_KpKm" ,100,0.9,2.2);
      book(_h8_Kpeta,"h08_Kpeta",100,1.3,2.7);
      book(_h8_Kmeta,"h08_Kmeta",100,1.3,2.7);
      book(_dalitz[7], "dalitz08",50,1.5,6.5,50,1.5,6.5);
      // KS0 KS0 eta'
      book(_h9_KS0KS0,"h09_KS0KS0",100,0.9,2.2);
      book(_h9_KS0eta,"h09_KS0eta",100,1.3,2.7);
      book(_dalitz[8], "dalitz09",50,1.5,6.5,50,1.5,6.5);
      // KL0 KL0 eta'
      book(_h10_KL0KL0,"h10_KL0KL0",100,0.9,2.2);
      book(_h10_KL0eta,"h10_KL0eta",100,1.3,2.7);
      book(_dalitz[9], "dalitz10",50,1.5,6.5,50,1.5,6.5);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip, Particles & pim, Particles & pi0,
			   Particles & Kp , Particles & Km ,
			   Particles & KS0, Particles & KL0,
			   Particles & eta, Particles & etaP) {
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
	else if (id == PID::K0S) {
          ++nstable;
	  KS0.push_back(p);
        }
	else if (id == PID::K0L) {
	  KL0.push_back(p);
          ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, pi0, Kp , Km, KS0, KL0, eta,etaP);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 441 )) {
	if(meson.mass()<2.93 || meson.mass()>3.03) continue;
	unsigned int nstable(0);
	Particles pip, pim, pi0, Kp , Km, KS0, KL0, eta,etaP;
	findDecayProducts(meson, nstable, pip, pim, pi0, Kp , Km, KS0, KL0, eta, etaP);
	if(nstable !=3) continue;
	// pi+ pi- eta
	if (pim.size()==1&&pip.size()==1&&eta.size()==1) {
	  double mplus  = (pip[0].momentum()+eta[0].momentum()).mass2();
	  double mminus = (pim[0].momentum()+eta[0].momentum()).mass2();
	  double mpipi  = (pip[0].momentum()+pim[0].momentum()).mass2();
	  _h1_pippim->fill(sqrt(mpipi));
	  _h1_pipeta->fill(sqrt(mplus));
	  _h1_pimeta->fill(sqrt(mminus));
	  _dalitz[0]->fill(mplus,mminus);
	}
	// pi0 pi0 eta
	else if (pi0.size()==2&&eta.size()==1) {
	  double m1    = (pi0[0].momentum()+eta[0].momentum()).mass2();
	  double m2    = (pi0[1].momentum()+eta[0].momentum()).mass2();
	  double mpipi = (pi0[0].momentum()+pi0[1].momentum()).mass2();
	  _h2_pi0pi0->fill(sqrt(mpipi));
	  _h2_pi0eta->fill(sqrt(m1));
	  _h2_pi0eta->fill(sqrt(m2));
	  _dalitz[1]->fill(m1,m2);
	  _dalitz[1]->fill(m2,m1);
	}
	// pi+ pi- eta
	else if (pim.size()==1&&pip.size()==1&&etaP.size()==1) {
	  double mplus  = (pip[0].momentum()+etaP[0].momentum()).mass2();
	  double mminus = (pim[0].momentum()+etaP[0].momentum()).mass2();
	  double mpipi  = (pip[0].momentum()+pim[0].momentum()).mass2();
	  _h3_pippim->fill(sqrt(mpipi));
	  _h3_pipeta->fill(sqrt(mplus));
	  _h3_pimeta->fill(sqrt(mminus));
	  _dalitz[2]->fill(mplus,mminus);
	}
	// pi0 pi0 eta'
	else if (pi0.size()==2&&etaP.size()==1) {
	  double m1    = (pi0[0].momentum()+etaP[0].momentum()).mass2();
	  double m2    = (pi0[1].momentum()+etaP[0].momentum()).mass2();
	  double mpipi = (pi0[0].momentum()+pi0[1].momentum()).mass2();
	  _h4_pi0pi0->fill(sqrt(mpipi));
	  _h4_pi0eta->fill(sqrt(m1));
	  _h4_pi0eta->fill(sqrt(m2));
	  _dalitz[3]->fill(m1,m2);
	  _dalitz[3]->fill(m2,m1);
	}
	// // K+K- eta
	else if (Km.size()==1&&Kp.size()==1&&eta.size()==1) {
	  double mplus  = (Kp[0].momentum()+eta[0].momentum()).mass2();
	  double mminus = (Km[0].momentum()+eta[0].momentum()).mass2();
	  double mKK    = (Kp[0].momentum()+Km [0].momentum()).mass2();
	  _h5_KpKm->fill(sqrt(mKK));
	  _h5_Kpeta->fill(sqrt(mplus));
	  _h5_Kmeta->fill(sqrt(mminus));
	  _dalitz[4]->fill(mplus,mminus);
	}
	// KS0 KS0 eta
	else if (KS0.size()==2&&eta.size()==1) {
	  double m1    = (KS0[0].momentum()+eta[0].momentum()).mass2();
	  double m2    = (KS0[1].momentum()+eta[0].momentum()).mass2();
	  double mKSKS = (KS0[0].momentum()+KS0[1].momentum()).mass2();
	  _h6_KS0KS0->fill(sqrt(mKSKS));
	  _h6_KS0eta->fill(sqrt(m1));
	  _h6_KS0eta->fill(sqrt(m2));
	  _dalitz[5]->fill(m1,m2);
	  _dalitz[5]->fill(m2,m1);
	}
	// KL0 KL0 eta
	else if (KL0.size()==2&&eta.size()==1) {
	  double m1    = (KL0[0].momentum()+eta[0].momentum()).mass2();
	  double m2    = (KL0[1].momentum()+eta[0].momentum()).mass2();
	  double mKLKL = (KL0[0].momentum()+KL0[1].momentum()).mass2();
	  _h7_KL0KL0->fill(sqrt(mKLKL));
	  _h7_KL0eta->fill(sqrt(m1));
	  _h7_KL0eta->fill(sqrt(m2));
	  _dalitz[6]->fill(m1,m2);
	  _dalitz[6]->fill(m2,m1);
	}
	// K+K- eta'
	else if (Km.size()==1&&Kp.size()==1&&etaP.size()==1) {
	  double mplus  = (Kp[0].momentum()+etaP[0].momentum()).mass2();
	  double mminus = (Km[0].momentum()+etaP[0].momentum()).mass2();
	  double mKK    = (Kp[0].momentum()+Km [0].momentum()).mass2();
	  _h8_KpKm->fill(sqrt(mKK));
	  _h8_Kpeta->fill(sqrt(mplus));
	  _h8_Kmeta->fill(sqrt(mminus));
	  _dalitz[7]->fill(mplus,mminus);
	}
	// KS0 KS0 eta'
	else if (KS0.size()==2&&etaP.size()==1) {
	  double m1    = (KS0[0].momentum()+etaP[0].momentum()).mass2();
	  double m2    = (KS0[1].momentum()+etaP[0].momentum()).mass2();
	  double mKSKS = (KS0[0].momentum()+KS0[1].momentum()).mass2();
	  _h9_KS0KS0->fill(sqrt(mKSKS));
	  _h9_KS0eta->fill(sqrt(m1));
	  _h9_KS0eta->fill(sqrt(m2));
	  _dalitz[8]->fill(m1,m2);
	  _dalitz[8]->fill(m2,m1);
	}
	// KL0 KL0 eta'
	else if (KL0.size()==2&&etaP.size()==1) {
	  double m1    = (KL0[0].momentum()+etaP[0].momentum()).mass2();
	  double m2    = (KL0[1].momentum()+etaP[0].momentum()).mass2();
	  double mKLKL = (KL0[0].momentum()+KL0[1].momentum()).mass2();
	  _h10_KL0KL0->fill(sqrt(mKLKL));
	  _h10_KL0eta->fill(sqrt(m1));
	  _h10_KL0eta->fill(sqrt(m2));
	  _dalitz[9]->fill(m1,m2);
	  _dalitz[9]->fill(m2,m1);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h1_pippim);
      normalize(_h1_pipeta);
      normalize(_h1_pimeta);
      normalize(_h2_pi0pi0);
      normalize(_h2_pi0eta);
      normalize(_h3_pippim);
      normalize(_h3_pipeta);
      normalize(_h3_pimeta);
      normalize(_h4_pi0pi0);
      normalize(_h4_pi0eta);
      normalize(_h5_KpKm);
      normalize(_h5_Kpeta);
      normalize(_h5_Kmeta);
      normalize(_h6_KS0KS0);
      normalize(_h6_KS0eta);
      normalize(_h7_KL0KL0);
      normalize(_h7_KL0eta);
      normalize(_h8_KpKm);
      normalize(_h8_Kpeta);
      normalize(_h8_Kmeta);
      normalize(_h9_KS0KS0);
      normalize(_h9_KS0eta);
      normalize(_h10_KL0KL0);
      normalize(_h10_KL0eta);
      for(unsigned int ix=0;ix<10;++ix)
	normalize(_dalitz[ix]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h1_pippim,_h1_pipeta,_h1_pimeta;
    Histo1DPtr _h2_pi0pi0,_h2_pi0eta;
    Histo1DPtr _h3_pippim,_h3_pipeta,_h3_pimeta;
    Histo1DPtr _h4_pi0pi0,_h4_pi0eta;
    Histo1DPtr _h5_KpKm,_h5_Kpeta,_h5_Kmeta;
    Histo1DPtr _h6_KS0KS0,_h6_KS0eta;
    Histo1DPtr _h7_KL0KL0,_h7_KL0eta;
    Histo1DPtr _h8_KpKm,_h8_Kpeta,_h8_Kmeta;
    Histo1DPtr _h9_KS0KS0,_h9_KS0eta;
    Histo1DPtr _h10_KL0KL0,_h10_KL0eta;
    
    Histo2DPtr _dalitz[10];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(MC_EtaC_Dalitz);

}
