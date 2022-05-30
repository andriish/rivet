// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  D+ -> pi+pi+pi-
  class E791_2001_I530320 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(E791_2001_I530320);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_h_pi,1,1,1);
      book(_dalitz, "dalitz",50,0.,3.1,50,0.0,3.1);
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
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 411 )) {
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
	  _h_pi->fill(m1);
	  _h_pi->fill(m2);
	  _dalitz->fill(m1,m2);
	  _dalitz->fill(m2,m1);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pi);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pi;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(E791_2001_I530320);

}
