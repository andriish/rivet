// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D_s+ -> KS0 KS0 pi+
  class BESIII_2022_I1945692 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I1945692);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_KK ,1,1,1);
      book(_h_Kpi,1,1,2);
      book(_dalitz,"dalitz",50,0.3,2.3,50,0.3,2.3);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , Particles & pi0,
			   Particles & K0) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::KPLUS || id == PID::KMINUS ||
	     id == PID::ETA   || id == PID::K0L) {
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
	  K0.push_back(p);
          ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, pi0, K0);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 431 )) {
      	unsigned int nstable(0);
      	Particles pip, pim, pi0, K0;
      	findDecayProducts(meson, nstable, pip, pim, pi0, K0);
      	if(nstable !=3) continue;
      	if(meson.pid()<0) {
      	  swap(pim,pip);
      	}
      	if (K0.size()==2&&pip.size()==1) {
     	  double m1  = (K0[0].momentum()+pip[0].momentum()).mass2();
     	  double m2  = (K0[1].momentum()+pip[0].momentum()).mass2();
     	  double mKK = (K0[0].momentum()+K0 [1].momentum()).mass2();
      	  _dalitz->fill(m1,m2);
      	  _dalitz->fill(m2,m1);
      	  _h_KK->fill(sqrt(mKK));
	  _h_Kpi->fill(sqrt(m1));
	  _h_Kpi->fill(sqrt(m2));
      	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Kpi);
      normalize(_h_KK );
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kpi,_h_KK;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I1945692);

}
