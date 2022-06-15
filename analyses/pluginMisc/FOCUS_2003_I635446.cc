// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D_s+/D+ -> pi+pi+pi-
  class FOCUS_2003_I635446 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(FOCUS_2003_I635446);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_h_pippim[0],1,1,1);
      book(_h_pippim[1],1,1,2);
      book(_h_pippim[2],2,1,1);
      book(_h_pippim[3],2,1,2);
      book(_dalitz[0], "dalitz1",50,0.,2.0,50,0.0,3.5);
      book(_dalitz[1], "dalitz2",50,0.,1.8,50,0.0,3.1);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::KPLUS || id == PID::KMINUS ||
	     id == PID::K0S   || id == PID::K0L ||
	     id == PID::PI0 ) {
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
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 411 ||
										   Cuts::abspid== 431 )) {
	unsigned int nstable(0);
	Particles pip, pim;
	findDecayProducts(meson, nstable, pip, pim);
	if(nstable !=3) continue;
	if(meson.pid()<0) {
	  swap(pim,pip);
	}
	if (pim.size()==1&&pip.size()==2) {
	  double m1 = (pim[0].momentum()+pip[0].momentum()).mass2();
	  double m2 = (pim[0].momentum()+pip[1].momentum()).mass2();
	  if(m1>m2) swap(m1,m2);
	  if(meson.abspid()==431) {
	    _dalitz[0]->fill(m1,m2);
	    _h_pippim[0]->fill(m1);
	    _h_pippim[1]->fill(m2);
	  }
	  else {
	    _dalitz[1]->fill(m1,m2);
	    _h_pippim[2]->fill(m1);
	    _h_pippim[3]->fill(m2);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_dalitz[ix]);
	normalize(_h_pippim[ix  ]);
	normalize(_h_pippim[ix+2]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pippim[4];
    Histo2DPtr _dalitz[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(FOCUS_2003_I635446);

}
