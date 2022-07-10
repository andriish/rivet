// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D0 -> eta eta pi0
  class BESIII_2018_I1662665 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1662665);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_pieta,1,1,1);
      book(_dalitz,"dalitz",50.,0.4,1.6,50,0.4,1.8);
    }

    
    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pi0, Particles & eta) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::KPLUS || id == PID::KMINUS ||
	     id == PID::K0S||id == PID::K0L ||
	     id == PID::PIPLUS || id == PID::PIMINUS) {
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
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pi0, eta);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 421 )) {
	unsigned int nstable(0);
	Particles pi0, eta;
	findDecayProducts(meson, nstable, pi0, eta);
	if(nstable !=3) continue;
	if (eta.size()==2&&pi0.size()==1) {
	  double m1 = (eta[0].momentum()+pi0[0].momentum()).mass2();
	  double m2 = (eta[1].momentum()+pi0[0].momentum()).mass2();
	  if(m1>m2) swap(m1,m2);
	  _dalitz->fill(m1,m2);
	  _h_pieta->fill(sqrt(m1));
	  _h_pieta->fill(sqrt(m2));
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pieta);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pieta;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2018_I1662665);

}
