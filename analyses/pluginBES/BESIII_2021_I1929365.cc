// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D_s+ -> pi+ pi0 pi0
  class BESIII_2021_I1929365 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1929365);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_h_pi0pi0,1,1,1);
      book(_h_pippi0,1,1,2);
      book(_dalitz, "dalitz",50,0.,3.5,50,0.,3.5);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pi0 , 
			   Particles & pip  , Particles & pim) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::PIPLUS ) {
       	  pip.push_back(p);
	  ++nstable;
	}
	else if (id == PID::PIMINUS ) {
	  pim.push_back(p);
	  ++nstable;
	}
	else if (id == PID::PI0) {
	  pi0.push_back(p);
	  ++nstable;
	}
	else if (id == PID::K0S || id == PID::KPLUS || id == PID::KMINUS) {
          ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pi0, pip , pim);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").
	    particles(Cuts::abspid== 431)) {
	unsigned int nstable(0);
	Particles pi0, pip , pim;
	findDecayProducts(meson, nstable, pi0, pip , pim);
	if(nstable !=3) continue;
	if(meson.pid()<0) {
	  swap(pip,pim);
	}
	if(pi0.size()==2&&pip.size()==1) {
	  double mp1  = (pip[0].momentum() +pi0[0].momentum()).mass2();
	  double mp2  = (pip[0].momentum() +pi0[1].momentum()).mass2();
	  double m0   = (pi0[0].momentum() +pi0[1].momentum()).mass();
	  _dalitz ->fill(mp1,mp2);
	  _dalitz ->fill(mp2,mp1);
	  // K_S0 veto
	  if(m0>0.458 && m0<0.520) continue;
	  _h_pippi0->fill(sqrt(mp1));
	  _h_pippi0->fill(sqrt(mp2));
	  _h_pi0pi0->fill(m0 );
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pi0pi0);
      normalize(_h_pippi0);
      normalize(_dalitz  );
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pi0pi0,_h_pippi0;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1929365);

}
