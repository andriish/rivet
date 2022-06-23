// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Ds+ -> K+ pi+ pi-
  class BESIII_2022_I2084294 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2084294);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_Kpip,1,1,1);
      book(_h_Kpim,1,1,2);
      book(_h_pipi,1,1,3);
      book(_dalitz,"dalitz",50,0.,2.3,50,0.3,3.5);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & Kp , Particles & Km ,
			   Particles & pip, Particles & pim) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::PIPLUS) {
	  pip.push_back(p);
	  ++nstable;
	}
	else if ( id == PID::PIMINUS) {
	  pim.push_back(p);
	  ++nstable;
	}
	else if (id == PID::ETA || id == PID::K0L ||
		 id == PID::PI0 || id == PID::K0S) {
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
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, Kp, Km, pip, pim);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 431 )) {
      	unsigned int nstable(0);
      	Particles Kp, Km, pip, pim;
      	findDecayProducts(meson, nstable, Kp, Km, pip, pim);
      	if(nstable !=3) continue;
      	if(meson.pid()<0) {
      	  swap(Km,Kp);
	  swap(pim,pip);
      	}
      	if (Kp.size()==1&&pip.size()==1&&pim.size()==1) {
     	  double mp    = (Kp [0].momentum()+pip[0].momentum()).mass2();
     	  double mm    = (Kp [0].momentum()+pim[0].momentum()).mass2();
     	  double mpipi = (pip[0].momentum()+pim[0].momentum()).mass2();
	  _dalitz->fill(mpipi,mm);
	  if(mpipi<sqr(.4676) || mpipi>sqr(0.5276)) {
	    _h_pipi->fill(sqrt(mpipi));
	    _h_Kpip->fill(sqrt(mp));
	    _h_Kpim->fill(sqrt(mm));
	  }
      	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Kpip);
      normalize(_h_Kpim);
      normalize(_h_pipi);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kpip,_h_Kpim,_h_pipi;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2084294);

}
