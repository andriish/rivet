// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BESIII_2017_I1621266 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2017_I1621266);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_h_pipi [0],1,1,1);
      book(_h_pieta[0],1,1,2);
      book(_h_pipi [1],1,1,3);
      book(_h_pieta[1],1,1,4);
      book(_dalitz[0], "dalitz_Jpsi" ,50,1., 9.,50,1., 9.);
      book(_dalitz[1], "dalitz_psi2S",50,1.,14.,50,1.,14.);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , Particles & eta) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::KPLUS || id == PID::KMINUS ||
	     id == PID::K0S) {
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
	else if (id == PID::ETAPRIME) {
	  eta.push_back(p);
	  ++nstable;
	}
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, eta);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").
	    particles(Cuts::pid== 443 || Cuts::pid==100443)) {
	unsigned int nstable(0);
	Particles pip, pim, eta;
	findDecayProducts(meson, nstable, pip, pim, eta);
	if(nstable !=3) continue;
	if(pip.size()==1 && pim.size()==1 && eta.size()==1) {
	  unsigned int iloc = meson.pid()/100000;
	  double mminus = (pim[0].momentum()+eta[0].momentum()).mass2();
	  double mplus  = (pip[0].momentum()+eta[0].momentum()).mass2();
	  double mneut  = (pip[0].momentum()+pim[0].momentum()).mass2();
	  _h_pipi[iloc] ->fill(sqrt(mneut ));
	  _h_pieta[iloc]->fill(sqrt(mplus ));
	  //_h_pieta[iloc]->fill(sqrt(mminus));
	  _dalitz[iloc]->fill(mplus,mminus);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_h_pipi [ix]);
	normalize(_h_pieta[ix]);
	normalize(_dalitz [ix]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pipi[2],_h_pieta[2];
    Histo2DPtr _dalitz[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2017_I1621266);

}
