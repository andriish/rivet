// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  D0 -> K- pi+ pi0
  class CLEOII_2001_I537154 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEOII_2001_I537154);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      
      book(_h_Kpi0,1,1,1);
      book(_h_Kpip,1,2,1);
      book(_h_pipi,1,3,1);

    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , Particles & pi0  ,
			   Particles & Kp  , Particles & Km  , Particles & K0) {
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
	  K0.push_back(p);
          ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, pi0, Kp , Km, K0);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const double E0   = 22.1e-5  ;
      static const double Ex   = -6.89e-5 ;
      static const double Ey   = -27.1e-5 ; 
      static const double Ex2  = 10.4e-5  ;
      static const double Exy  = 38.2e-5  ;
      static const double Ey2  = 12.4e-5  ;
      static const double Ex3  = -3.00e-5;
      static const double Ex2y = -7.97e-5;
      static const double Exy2 = -12.8e-5;
      static const double Ey3  = -0.53e-5;
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 421)) {
	unsigned int nstable(0);
	Particles pip, pim, pi0, Kp , Km, K0;
	findDecayProducts(meson, nstable, pip, pim, pi0, Kp , Km, K0);
	if(meson.pid()<0) {
	  swap(pim,pip);
	  swap(Kp,Km);
	}
	if(pip.size()==1&&Km.size()==1&&pi0.size()==1) {
	  double mKpip = (Km [0].momentum() + pip[0].momentum()).mass2();
	  double mKpi0 = (Km [0].momentum() + pi0[0].momentum()).mass2();
	  double mpipi = (pi0[0].momentum() + pip[0].momentum()).mass2();
	  double eff = E0 + Ex*mKpip + Ey*mpipi +Ex2*sqr(mKpip)+Exy*mKpip*mpipi+Ey2*sqr(mpipi)
	    +Ex3*pow(mKpip,3)+Ex2y*sqr(mKpip)*mpipi +Exy2*mKpip*sqr(mpipi)+Ey3*pow(mpipi,3);
	  _h_Kpi0->fill(mKpi0,eff);
	  _h_Kpip->fill(mKpip,eff);
	  _h_pipi->fill(mpipi,eff);
	}
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Kpi0);
      normalize(_h_Kpip);
      normalize(_h_pipi);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_Kpi0,_h_Kpip,_h_pipi;
    //@}


  };


  DECLARE_RIVET_PLUGIN(CLEOII_2001_I537154);

}
