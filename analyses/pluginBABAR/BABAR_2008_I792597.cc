// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D_s+ -> pi+pi+pi-
  class BABAR_2008_I792597 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2008_I792597);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_h_pippim[0],1,1,3);
      book(_h_pippip   ,1,1,4);
      book(_h_pippim[1],1,1,1);
      book(_h_pippim[2],1,1,2);
      book(_dalitz, "dalitz",50,0.,3.5,50,0.0,3.5);
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
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 431 )) {
	unsigned int nstable(0);
	Particles pip, pim;
	findDecayProducts(meson, nstable, pip, pim);
	if(nstable !=3) continue;
	if(meson.pid()<0) {
	  swap(pim,pip);
	}
	if (pim.size()==1&&pip.size()==2) {
	  // kinematic variables
	  double x[3] = {(pim[0].momentum()+pip[0].momentum()).mass2(),
			 (pim[0].momentum()+pip[1].momentum()).mass2(),
			 (pip[0].momentum()+pip[1].momentum()).mass2()};
	  if(x[0]>x[1]) swap(x[0],x[1]);
	  _dalitz->fill(x[0],x[1]);
	  _dalitz->fill(x[1],x[0]);
	  // fill plots
	  _h_pippim[0]->fill(x[0]);
	  _h_pippim[0]->fill(x[1]);
	  _h_pippim[1]->fill(x[0]);
	  _h_pippim[2]->fill(x[1]);
	  _h_pippip   ->fill(x[2]);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	normalize(_h_pippim[ix]);
      normalize(_h_pippip);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pippim[3],_h_pippip;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2008_I792597);

}
