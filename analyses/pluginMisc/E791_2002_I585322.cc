// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D+ -> K-pi+pi+
  class E791_2002_I585322 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(E791_2002_I585322);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_h_low,1,1,1);
      book(_h_high,1,1,2);
      book(_dalitz, "dalitz",50,0.,3.1,50,0.3,3.1);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim,
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
	else if (id == PID::PIPLUS) {
	  pip.push_back(p);
	  ++nstable;
	}
	else if (id == PID::PIMINUS) {
	  pim.push_back(p);
	  ++nstable;
	}
	else if (id == PID::PI0 || id == PID::K0S) {
          ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, Kp , Km);
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
	Particles pip, pim, Kp , Km;
	findDecayProducts(meson, nstable, pip, pim, Kp , Km);
	if(nstable !=3) continue;
	if(meson.pid()<0) {
	  swap(pim,pip);
	  swap(Kp,Km);
	}
	if(pip.size()==2&&Km.size()==1) {
	  double mplus  = (Km[0].momentum() +pip[0].momentum()).mass2();
	  double mminus = (Km[0].momentum() +pip[1].momentum()).mass2();
	  if(mplus<mminus) swap(mplus,mminus);
	  _h_low ->fill(mminus);
	  _h_high->fill(mplus );
	  _dalitz->fill(mminus,mplus );
	  _dalitz->fill(mplus ,mminus);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_low );
      normalize(_h_high);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_high,_h_low;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(E791_2002_I585322);

}
