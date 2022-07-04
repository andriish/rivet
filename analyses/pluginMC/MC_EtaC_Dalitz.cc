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
      book(_dalitz[4], "dalitz05",50,1.,7.,50,1.,7.);
      // KS0 KS0 eta
      book(_h6_KS0KS0,"h06_KS0KS0",100,0.5,2.5);
      book(_h6_KS0eta,"h06_KS0eta",100,1.0,3.0);
      book(_dalitz[5], "dalitz06",50,1.,7.,50,1.,7.);
      // KL0 KL0 eta
      book(_h7_KL0KL0,"h07_KL0KL0",100,0.5,2.5);
      book(_h7_KL0eta,"h07_KL0eta",100,1.0,3.0);
      book(_dalitz[6], "dalitz07",50,1.,7.,50,1.,7.);
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
      book(_dalitz[9] ,  "dalitz10",50,1.5,6.5,50,1.5,6.5);
      // K+K- pi0
      book(_h11_KpKm ,"h11_KpKm" ,100,0.9,2.8);
      book(_h11_Kppi0,"h11_Kppi0",100,0.6,2.8);
      book(_h11_Kmpi0,"h11_Kmpi0",100,0.6,2.8);
      book(_dalitz[10],"dalitz11",50,0.3,6.5,50,0.3,6.5);
      // KS0 KS0 pi0
      book(_h12_KS0KS0,"h12_KS0KS0",100,0.9,2.8);
      book(_h12_KS0pi0,"h12_KS0pi0",100,0.6,2.8);
      book(_dalitz[11] ,"dalitz12",50,0.3,6.5,50,0.3,6.5);
      // KL0 KL0 pi0
      book(_h13_KL0KL0,"h13_KL0KL0",100,0.9,2.8);
      book(_h13_KL0pi0,"h13_KL0pi0",100,0.6,2.8);
      book(_dalitz[12] ,  "dalitz13",50,0.3,6.5,50,0.3,6.5);
      // KS0 K+ pi-
      book(_h14_KpKS0 ,"h14_KpKS0" ,100,0.9,2.8);
      book(_h14_Kppim ,"h14_Kppim" ,100,0.6,2.8);
      book(_h14_KS0pim,"h14_KS0pim",100,0.6,2.8);
      book(_dalitz[13] ,"dalitz14",50,0.3,6.5,50,0.3,6.5);
      // KS0 K- pi+
      book(_h15_KmKS0 ,"h15_KmKS0" ,100,0.9,2.8);
      book(_h15_Kmpip ,"h15_Kmpip" ,100,0.6,2.8);
      book(_h15_KS0pip,"h15_KS0pip",100,0.6,2.8);
      book(_dalitz[14] ,  "dalitz15",50,0.3,6.5,50,0.3,6.5);
      // KL0 K+ pi-
      book(_h16_KpKL0 ,"h16_KpKL0" ,100,0.9,2.8);
      book(_h16_Kppim ,"h16_Kppim" ,100,0.6,2.8);
      book(_h16_KL0pim,"h16_KL0pim",100,0.6,2.8);
      book(_dalitz[15] ,"dalitz16",50,0.3,6.5,50,0.3,6.5);
      // KL0 K- pi+
      book(_h17_KmKL0 ,"h17_KmKL0" ,100,0.9,2.8);
      book(_h17_Kmpip ,"h17_Kmpip" ,100,0.6,2.8);
      book(_h17_KL0pip,"h17_KL0pip",100,0.6,2.8);
      book(_dalitz[16] ,"dalitz17",50,0.3,6.5,50,0.3,6.5);
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
	// K+ K- pi0
	else if (Km.size()==1&&Kp.size()==1&&pi0.size()==1) {
	  double mplus  = (Kp[0].momentum()+pi0[0].momentum()).mass2();
	  double mminus = (Km[0].momentum()+pi0[0].momentum()).mass2();
	  double mKK    = (Kp[0].momentum()+Km [0].momentum()).mass2();
	  _h11_KpKm->fill(sqrt(mKK));
	  _h11_Kppi0->fill(sqrt(mplus));
	  _h11_Kmpi0->fill(sqrt(mminus));
	  _dalitz[10]->fill(mplus,mminus);
	}
	// KS0 KS0 pi0
	else if (KS0.size()==2&&pi0.size()==1) {
	  double mplus  = (KS0[0].momentum()+pi0[0].momentum()).mass2();
	  double mminus = (KS0[1].momentum()+pi0[0].momentum()).mass2();
	  double mKK    = (KS0[0].momentum()+KS0[1].momentum()).mass2();
	  _h12_KS0KS0->fill(sqrt(mKK));
	  _h12_KS0pi0->fill(sqrt(mplus));
	  _h12_KS0pi0->fill(sqrt(mminus));
	  _dalitz[11]->fill(mplus,mminus);
	}
	// KL0 KL0 pi0
	else if (KL0.size()==2&&pi0.size()==1) {
	  double mplus  = (KL0[0].momentum()+pi0[0].momentum()).mass2();
	  double mminus = (KL0[1].momentum()+pi0[0].momentum()).mass2();
	  double mKK    = (KL0[0].momentum()+KL0[1].momentum()).mass2();
	  _h13_KL0KL0->fill(sqrt(mKK));
	  _h13_KL0pi0->fill(sqrt(mplus));
	  _h13_KL0pi0->fill(sqrt(mminus));
	  _dalitz[12]->fill(mplus,mminus);
	}
	// KS0 K+ pi-
	else if (KS0.size()==1&&Kp.size()==1&&pim.size()==1) {
	  double mplus  = (Kp[0].momentum()+pim[0].momentum()).mass2();
	  double mminus = (KS0[0].momentum()+pim[0].momentum()).mass2();
	  double mKK    = (KS0[0].momentum()+Kp[0].momentum()).mass2();
	  _h14_KpKS0->fill(sqrt(mKK));
	  _h14_Kppim->fill(sqrt(mplus));
	  _h14_KS0pim->fill(sqrt(mminus));
	  _dalitz[13]->fill(mplus,mminus);
	}
	// KS0 K- pip
	else if (KS0.size()==1&&Km.size()==1&&pip.size()==1) {
	  double mplus  = (Km[0].momentum()+pip[0].momentum()).mass2();
	  double mminus = (KS0[0].momentum()+pip[0].momentum()).mass2();
	  double mKK    = (KS0[0].momentum()+Km[0].momentum()).mass2();
	  _h15_KmKS0->fill(sqrt(mKK));
	  _h15_Kmpip->fill(sqrt(mplus));
	  _h15_KS0pip->fill(sqrt(mminus));
	  _dalitz[14]->fill(mplus,mminus);
	}
	// KL0 K+ pi-
	else if (KL0.size()==1&&Kp.size()==1&&pim.size()==1) {
	  double mplus  = (Kp[0].momentum()+pim[0].momentum()).mass2();
	  double mminus = (KL0[0].momentum()+pim[0].momentum()).mass2();
	  double mKK    = (KL0[0].momentum()+Kp[0].momentum()).mass2();
	  _h16_KpKL0->fill(sqrt(mKK));
	  _h16_Kppim->fill(sqrt(mplus));
	  _h16_KL0pim->fill(sqrt(mminus));
	  _dalitz[15]->fill(mplus,mminus);
	}
	// KL0 K- pip
	else if (KL0.size()==1&&Km.size()==1&&pip.size()==1) {
	  double mplus  = (Km[0].momentum()+pip[0].momentum()).mass2();
	  double mminus = (KL0[0].momentum()+pip[0].momentum()).mass2();
	  double mKK    = (KL0[0].momentum()+Km[0].momentum()).mass2();
	  _h17_KmKL0->fill(sqrt(mKK));
	  _h17_Kmpip->fill(sqrt(mplus));
	  _h17_KL0pip->fill(sqrt(mminus));
	  _dalitz[16]->fill(mplus,mminus);
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
      normalize(_h11_KpKm );
      normalize(_h11_Kppi0);
      normalize(_h11_Kmpi0);
      normalize(_h12_KS0KS0);
      normalize(_h12_KS0pi0);
      normalize(_h13_KL0KL0);
      normalize(_h13_KL0pi0);
      normalize(_h14_KpKS0);
      normalize(_h14_Kppim);
      normalize(_h14_KS0pim);
      normalize(_h15_KmKS0);
      normalize(_h15_Kmpip);
      normalize(_h15_KS0pip);
      normalize(_h16_KpKL0);
      normalize(_h16_Kppim);
      normalize(_h16_KL0pim);
      normalize(_h17_KmKL0);
      normalize(_h17_Kmpip);
      normalize(_h17_KL0pip);
      for(unsigned int ix=0;ix<17;++ix)
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
    Histo1DPtr _h11_KpKm,_h11_Kppi0,_h11_Kmpi0;
    Histo1DPtr _h12_KS0KS0,_h12_KS0pi0;
    Histo1DPtr _h13_KL0KL0,_h13_KL0pi0;
    Histo1DPtr _h14_KpKS0,_h14_Kppim,_h14_KS0pim;
    Histo1DPtr _h15_KmKS0,_h15_Kmpip,_h15_KS0pip;
    Histo1DPtr _h16_KpKL0,_h16_Kppim,_h16_KL0pim;
    Histo1DPtr _h17_KmKL0,_h17_Kmpip,_h17_KL0pip;
    Histo2DPtr _dalitz[17];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(MC_EtaC_Dalitz);

}
