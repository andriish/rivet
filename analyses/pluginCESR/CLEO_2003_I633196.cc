// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D0 -> KS0 pi+ pi-
  class CLEO_2003_I633196 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2003_I633196);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Histograms
      book(_h_Kpim,1,1,1);
      book(_h_pipi,1,1,2);
      book(_h_Kpip,1,1,3);
      book(_dalitz, "dalitz",50,0.,3.,50,0.,2.);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , Particles & K0) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::KPLUS || id == PID::KMINUS || id == PID::PI0 || id == PID::K0L) {
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
	else if (id == PID::K0S) {
	  K0.push_back(p);
          ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, K0);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").
	    particles(Cuts::abspid== 421)) {
	unsigned int nstable(0);
	Particles pip, pim, K0;
	findDecayProducts(meson, nstable, pip, pim, K0);
	if(nstable !=3) continue;
	if(meson.pid()<0) {
	  swap(pim,pip);
	}
	if(pim.size()==1&&pip.size()==1&&K0.size()==1) {
	  double mminus = (pim[0].momentum()+K0[0].momentum() ).mass2();
	  double mplus  = (pip[0].momentum()+K0[0].momentum() ).mass2();
	  double mpipi  = (pip[0].momentum()+pim[0].momentum()).mass2();
	  _h_Kpip->fill(mplus);
	  _h_Kpim->fill(mminus);
	  _h_pipi->fill(mpipi);
	  _dalitz->fill(mminus,mpipi); 
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Kpip);
      normalize(_h_pipi);
      normalize(_h_Kpim);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kpim,_h_pipi,_h_Kpip;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_2003_I633196);

}
