// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief psi(2S) -> J/psi pi+pi-
  class MARKII_1979_I144382 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MARKII_1979_I144382);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(),"UFS");
      book(_hist, 1, 1, 1);

    }


    void findDecayProducts(const Particle & mother,
			   unsigned int & nstable,
			   Particles& pip, Particles& pim,
			   Particles& pi0, Particles & onium) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
      	if ( id == PID::PIMINUS) {
	  pim.push_back(p);
	  ++nstable;
	}
       	else if (id == PID::PIPLUS) {
       	  pip.push_back(p);
       	  ++nstable;
       	}
       	else if (id == PID::PI0) {
       	  pi0.push_back(p);
       	  ++nstable;
       	}
	else if (abs(id)%1000==443 || abs(id)%1000==553) {
	  onium.push_back(p);
	  ++nstable;
	}
	else if ( !p.children().empty() ) {
	  findDecayProducts(p,nstable,pip,pim,pi0,onium);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over unstable particles
      for(const Particle& ups : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==100443)) {
      	unsigned int nstable(0);
      	Particles pip, pim, pi0, onium;
      	findDecayProducts(ups,nstable,pip,pim,pi0,onium);
      	// check for onium
      	if(onium.size() !=1 || onium[0].pid()!=443 || nstable !=3) continue;
      	// check for pipi
      	if( ! (pip.size()==1 && pim.size() ==1) ) continue;
      	FourMomentum q = pip[0].momentum()+pim[0].momentum();
      	_hist->fill(q.mass2());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_hist);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _hist;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(MARKII_1979_I144382);

}
