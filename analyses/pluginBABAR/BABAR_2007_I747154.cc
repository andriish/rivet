// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D0 -> pi+pi-pi0
  class BABAR_2007_I747154 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2007_I747154);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_h_pi[0],1,1,1);
      book(_h_pi[1],1,1,2);
      book(_dalitz, "dalitz",50,0.,3.2,50,0.0,3.2);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , Particles & pi0) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::KPLUS ) {
	  ++nstable;
	}
	else if (id == PID::KMINUS ) {
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
          ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, pi0);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 421 )) {
	unsigned int nstable(0);
	Particles pip, pim, pi0;
	findDecayProducts(meson, nstable, pip, pim, pi0);
	if(nstable !=3) continue;
	if(meson.pid()<0) {
	  swap(pim,pip);
	}
	if (nstable==3 && pim.size()==1&&pi0.size()==1&&pip.size()==1) {
	  double mplus  = (pip[0].momentum()+pi0[0].momentum()).mass2();
	  double mminus = (pim[0].momentum()+pi0[0].momentum()).mass2();
	  _h_pi[0]->fill(mplus);
	  _h_pi[1]->fill(mminus);
	  _dalitz->fill(mminus,mplus);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pi[0]);
      normalize(_h_pi[1]);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pi[2];
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2007_I747154);

}
