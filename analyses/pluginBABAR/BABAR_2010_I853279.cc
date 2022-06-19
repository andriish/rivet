// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D0 -> KS0 pi+pi- and K+K-
  class BABAR_2010_I853279 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2010_I853279);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // Histograms
      book(_h_Kpim,1,1,1);
      book(_h_Kpip,1,1,2);
      book(_h_pipi,1,1,3);
      book(_dalitz[0], "dalitz1",50,0.3,3.2,50,0.3,3.2);
      book(_h_K0Km,1,1,4);
      book(_h_K0Kp,1,1,5);
      book(_h_KpKm,1,1,6);
      book(_dalitz[1], "dalitz2",50,0.9,1.9,50,0.9,1.9);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , Particles & K0, Particles & Kp,  Particles & Km) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::KPLUS || id == PID::KMINUS || id == PID::PI0 || id == PID::K0L) {
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
	if(nstable !=3) continue;
	if(meson.pid()<0) {
	  swap(pim,pip);
	  swap(Km,Kp);
	}
	if(pim.size()==1&&pip.size()==1&&K0.size()==1) {
	  double mminus = (pim[0].momentum()+K0[0].momentum() ).mass2();
	  double mplus  = (pip[0].momentum()+K0[0].momentum() ).mass2();
	  double mpipi  = (pip[0].momentum()+pim[0].momentum()).mass2();
	  _h_Kpip->fill(mplus);
	  _h_Kpim->fill(mminus);
	  _h_pipi->fill(mpipi);
	  _dalitz[0]->fill(mminus,mplus); 
	}
	else if(Km.size()==1&&Kp.size()==1&&K0.size()==1) {
	  double mminus = (Km[0].momentum()+K0[0].momentum() ).mass2();
	  double mplus  = (Kp[0].momentum()+K0[0].momentum() ).mass2();
	  double mKK  = (Kp[0].momentum()+Km[0].momentum()).mass2();
	  _h_K0Kp->fill(mplus);
	  _h_K0Km->fill(mminus);
	  _h_KpKm->fill(mKK);
	  _dalitz[1]->fill(mplus,mKK); 
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Kpim);
      normalize(_h_Kpip);
      normalize(_h_pipi);
      normalize(_dalitz[0]);
      normalize(_h_K0Km);
      normalize(_h_K0Kp);
      normalize(_h_KpKm);
      normalize(_dalitz[1]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kpim,_h_pipi,_h_Kpip;
    Histo1DPtr _h_K0Km,_h_KpKm,_h_K0Kp;
    Histo2DPtr _dalitz[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2010_I853279);

}
