// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D+/Ds+ -> K+ pi+ pi-
  class FOCUS_2004_I654030 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(FOCUS_2004_I654030);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // histograms
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_Kpim[ix],1+ix,1,1);
	book(_h_pipi[ix],1+ix,1,2);
	book(_dalitz[ix],"dalitz_"+to_str(ix+1),50,0.3,3.5,50,0.,2.3);
      }
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
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 411 or
										   Cuts::abspid== 431 )) {
      	unsigned int nstable(0);
      	Particles Kp, Km, pip, pim;
      	findDecayProducts(meson, nstable, Kp, Km, pip, pim);
      	if(nstable !=3) continue;
      	if(meson.pid()<0) {
      	  swap(Km,Kp);
	  swap(pim,pip);
      	}
      	if (Kp.size()==1&&pip.size()==1&&pim.size()==1) {
     	  double mm    = (Kp [0].momentum()+pim[0].momentum()).mass2();
     	  double mpipi = (pip[0].momentum()+pim[0].momentum()).mass2();
	  unsigned int iloc = meson.abspid()==411 ? 0 : 1;
	  _dalitz[iloc]->fill(mm,mpipi);
	  _h_pipi[iloc]->fill(mpipi);
	  _h_Kpim[iloc]->fill(mm);
      	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_h_Kpim[ix]);
	normalize(_h_pipi[ix]);
	normalize(_dalitz[ix]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kpim[2],_h_pipi[2];
    Histo2DPtr _dalitz[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(FOCUS_2004_I654030);

}
