// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/psi and psi(2S) -> pi+pi-pi0
  class BESIII_2012_I1088606 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2012_I1088606);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_h_pipi[0],1,1,1);
      book(_h_pipi[1],1,1,2);
      book(_dalitz[0], "dalitz_Jpsi" ,50,0.,9.,50,0.0,9.);
      book(_dalitz[1], "dalitz_psi2S",50,0.,9.,50,0.0,9.);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , Particles & pi0) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::KPLUS || id == PID::KMINUS || id == PID::K0S || id == PID::K0L ) {
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
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, pi0);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").
	    particles(Cuts::pid== 443 or Cuts::pid==100443)) {
	unsigned int nstable(0);
	Particles pip, pim, pi0;
	findDecayProducts(meson, nstable, pip, pim, pi0);
	if(nstable !=3) continue;
	if(pip.size()==1 && pim.size()==1 && pi0.size()==1) {
	    double mminus = (pim[0].momentum()+pi0[0].momentum()).mass2();
	    double mplus  = (pip[0].momentum()+pi0[0].momentum()).mass2();
	    double mneut  = (pip[0].momentum()+pim[0].momentum()).mass2();
	    unsigned int iloc = meson.pid()/100000;
	    _h_pipi[iloc]->fill(sqrt(mneut ));
	    _h_pipi[iloc]->fill(sqrt(mplus ));
	    _h_pipi[iloc]->fill(sqrt(mminus));
	    _dalitz[iloc]->fill(mplus,mminus);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_h_pipi  [ix]);
	normalize(_dalitz[ix]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pipi[2];
    Histo2DPtr _dalitz[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2012_I1088606);

}
