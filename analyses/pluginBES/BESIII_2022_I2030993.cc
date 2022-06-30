// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D_s+ -> pi+ pi0 eta'
  class BESIII_2022_I2030993 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2030993);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_etapip,1,1,1);
      book(_h_etapi0,1,1,2);
      book(_h_pipi , 1,1,3);
      book(_dalitz, "dalitz",50,0.,1.1,50,1,3.5);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , Particles & pi0,
			   Particles & eta) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::KPLUS || id == PID::KMINUS ||
	     id == PID::K0S||id == PID::K0L) {
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
	else if (id == PID::ETAPRIME) {
	  eta.push_back(p);
          ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, pi0, eta);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 431 )) {
	unsigned int nstable(0);
	Particles pip, pim, pi0, eta;
	findDecayProducts(meson, nstable, pip, pim, pi0, eta);
	if(nstable !=3) continue;
	if(meson.pid()<0) {
	  swap(pim,pip);
	}
	if (eta.size()==1&&pi0.size()==1&&pip.size()==1) {
	  double mplus = (eta[0].momentum()+pip[0].momentum()).mass2();
	  double mzero = (eta[0].momentum()+pi0[0].momentum()).mass2();
	  double mpipi = (pip[0].momentum()+pi0[0].momentum()).mass2();
	  _dalitz  ->fill(mpipi,mplus);
	  _h_pipi  ->fill(sqrt(mpipi));
	  _h_etapip->fill(sqrt(mplus));
	  _h_etapi0->fill(sqrt(mzero));
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pipi  );
      normalize(_h_etapip);
      normalize(_h_etapi0);
      normalize(_dalitz  );
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pipi, _h_etapip,_h_etapi0;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2030993);

}
