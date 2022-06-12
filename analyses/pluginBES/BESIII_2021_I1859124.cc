// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D+ -> K+ KS0 pi0
  class BESIII_2021_I1859124 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1859124);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_h_Kppi0,1,1,1);
      book(_h_K0pi0,1,1,2);
      book(_h_KpK0,1,1,3);
      book(_dalitz, "dalitz",50,0.3,2.,50,0.3,2.);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pi0 , Particles & K0,
			   Particles & Kp  , Particles & Km) {
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
	else if (id == PID::K0S) {
	  K0.push_back(p);
	  ++nstable;
	}
	else if (id == PID::PI0) {
	  pi0.push_back(p);
	  ++nstable;
	}
	else if (id == PID::PIPLUS || id == PID::PIMINUS) {
          ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pi0, K0, Kp , Km);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").
	    particles(Cuts::abspid== 411)) {
	unsigned int nstable(0);
	Particles pi0, K0, Kp , Km;
	findDecayProducts(meson, nstable, pi0, K0, Kp , Km);
	if(nstable !=3) continue;
	if(meson.pid()<0) {
	  swap(Kp,Km);
	}
	if(pi0.size()==1&&Kp.size()==1&&K0.size()==1) {
	  double mplus  = (Kp[0].momentum() +pi0[0].momentum()).mass2();
	  double mneut  = (K0[0].momentum() +pi0[0].momentum()).mass2();
	  double mKK    = (Kp[0].momentum() + K0[0].momentum()).mass2();
	  _h_Kppi0->fill(sqrt(mplus));
	  _h_K0pi0->fill(sqrt(mneut));
	  _h_KpK0 ->fill(sqrt(mKK  ));
	  _dalitz ->fill(mplus,mneut);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Kppi0);
      normalize(_h_K0pi0);
      normalize(_h_KpK0 );
      normalize(_dalitz );
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kppi0,_h_K0pi0,_h_KpK0;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1859124);

}
