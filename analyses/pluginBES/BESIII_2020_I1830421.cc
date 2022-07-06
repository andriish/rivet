// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief eta' -> pi+pi- e+e- decay
  class BESIII_2020_I1830421 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2020_I1830421);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      //book(_h_mee  , 1, 1, 1);
      book(_h_mpipi, 1, 1, 2);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , 
			   Particles &  ep , Particles &  em) {
      for(const Particle & p : mother.children()) {
	int id = p.pid();
        if (id == PID::PIPLUS) {
	  pip.push_back(p);
	  ++nstable;
	}
	else if (id == PID::PIMINUS) {
	  pim.push_back(p);
	  ++nstable;
	}
	else if (id == PID::EMINUS) {
	  em.push_back(p);
	  ++nstable;
	}
	else if (id == PID::EPLUS) {
	  ep.push_back(p);
	  ++nstable;
	}
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, ep,em);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Loop over eta' mesons
      for (const Particle& p :  apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==331)) {
	unsigned int nstable=0;
	Particles pip, pim, ep, em;
	findDecayProducts(p,nstable,pip, pim, ep, em);
	if( nstable==4 &&
	    pip.size() == 1 && pim.size() == 1 &&
	    em.size()  == 1 && ep.size()  == 1 ) {
	  //_h_mee  ->fill((em [0].momentum()+ep [0].momentum()).mass());
	  _h_mpipi->fill((pim[0].momentum()+pip[0].momentum()).mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      //normalize(_h_mee);
      normalize(_h_mpipi);
    }

    /// @}


    /// @name Histograms
    /// @{
    //Histo1DPtr _h_mee;
    Histo1DPtr _h_mpipi;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2020_I1830421);

}
