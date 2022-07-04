// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/psi -> K+K-pi0
  class BESIII_2019_I1731057 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2019_I1731057);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_KpKm [ix],1+ix,1,1);
	book(_h_Kppi0[ix],1+ix,1,2);
      }
      book(_dalitz, "dalitz",50,0.,7.,50,0.0,7.);
    }


    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pi0, Particles & Kp, Particles & Km) {
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
	else if (id == PID::PIPLUS || id == PID::PIMINUS ||
		 id == PID::K0S    || id == PID::K0L) {
	  ++nstable;
	}
	else if (id == PID::PI0) {
	  pi0.push_back(p);
	  ++nstable;
	}
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pi0, Kp , Km);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").
	    particles(Cuts::pid== 443)) {
	unsigned int nstable(0);
	Particles pi0, Kp , Km;
	findDecayProducts(meson, nstable, pi0, Kp , Km);
	if(nstable !=3) continue;
	if(Kp.size()==1 && Km.size()==1 && pi0.size()==1) {
	    double mminus = (Km[0].momentum()+pi0[0].momentum()).mass2();
	    double mplus  = (Kp[0].momentum()+pi0[0].momentum()).mass2();
	    double mneut  = (Kp[0].momentum()+Km[0].momentum()).mass2();
	    _dalitz->fill(mplus,mminus);
	    mplus=sqrt(mplus);
	    mminus=sqrt(mminus);
	    mneut=sqrt(mneut);
	    _h_KpKm [0]->fill(mneut );
	    _h_Kppi0[0]->fill(mplus );
	    _h_Kppi0[0]->fill(mminus);
	    if(mplus>1.05 && mminus>1.05) {
	      _h_KpKm [1]->fill(mneut );
	      _h_Kppi0[1]->fill(mplus );
	      _h_Kppi0[1]->fill(mminus);
	    }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_h_KpKm [ix],1.,false);
	normalize(_h_Kppi0[ix],1.,false);
      }
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_KpKm[2],_h_Kppi0[2];
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2019_I1731057);

}
