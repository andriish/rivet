// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"


namespace Rivet {


  /// @brief J/psi dalitz decays
  class BABAR_2017_I1512302 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2017_I1512302);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_h_pippim,1,1,1);
      book(_h_pippi0,1,1,2);
      book(_dalitz_3pi, "dalitz_3pi",50,0.,9.,50,0.0,9.);
      book(_h_KpKm ,2,1,1);
      book(_h_Kppi0,2,1,2);
      book(_dalitz_KpKmpi, "dalitz_KpKmpi",50,0.,7.,50,0.0,7.);
      book(_h_K0Kp ,3,1,1);
      book(_h_K0pip,3,1,2);
      book(_h_Kppip,3,1,3);
      book(_dalitz_K0Kppim, "dalitz_K0Kppim",50,0.,8.,50,0.,8.);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , Particles & pi0  ,
			   Particles & Kp  , Particles & Km  , Particles & K0) {
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
	else if (id == PID::PI0) {
	  pi0.push_back(p);
	  ++nstable;
	}
	else if (id == PID::K0S) {
	  K0.push_back(p);
          ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, pi0, Kp , Km, K0);
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
	Particles pip, pim, pi0, Kp , Km, K0;
	findDecayProducts(meson, nstable, pip, pim, pi0, Kp , Km, K0);
	if(nstable !=3) continue;
	if(pip.size()==1 && pim.size()==1 && pi0.size()==1) {
	    double mminus = (pim[0].momentum()+pi0[0].momentum()).mass2();
	    double mplus  = (pip[0].momentum()+pi0[0].momentum()).mass2();
	    double mneut  = (pip[0].momentum()+pim[0].momentum()).mass2();
	    _h_pippim->fill(mneut );
	    _h_pippi0->fill(mplus );
	    _h_pippi0->fill(mminus);
	    _dalitz_3pi->fill(mplus,mminus);
	}
	else if(Kp.size()==1 && Km.size()==1 && pi0.size()==1) {
	    double mminus = (Km[0].momentum()+pi0[0].momentum()).mass2();
	    double mplus  = (Kp[0].momentum()+pi0[0].momentum()).mass2();
	    double mneut  = (Kp[0].momentum()+Km[0].momentum()).mass2();
	    _h_KpKm->fill(mneut );
	    _h_Kppi0->fill(mplus );
	    _h_Kppi0->fill(mminus);
	    _dalitz_KpKmpi->fill(mplus,mminus);
	}
	else if (((Km.size()==1&&pip.size()==1) ||
		  (Kp.size()==1&&pim.size()==1) )&&K0.size()==1&&
		 meson.mass()>2.922 && meson.mass()<3.039) {
	  if(Km.size()==1) {
	    swap(Km,Kp);
	    swap(pim,pip);
	  }
	  double mplus  = (Kp[0].momentum()  +  K0[0].momentum()).mass2();
	  double mminus = (K0[0].momentum()  + pim[0].momentum()).mass2();
	  double mKK    = (Kp[0].momentum()  +  K0[0].momentum()).mass2();
	  _h_K0Kp ->fill(mKK);
	  _h_Kppip->fill(mplus);
	  _h_K0pip->fill(mminus);
	  _dalitz_K0Kppim->fill(mplus,mminus);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pippim,1.,false);
      normalize(_h_pippi0,1.,false);
      normalize(_dalitz_3pi);
      normalize(_h_KpKm,1.,false);
      normalize(_h_Kppi0,1.,false);
      normalize(_dalitz_KpKmpi);
      normalize(_h_Kppip);
      normalize(_h_K0pip);
      normalize(_h_K0Kp );
      normalize(_dalitz_K0Kppim);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pippim,_h_pippi0;
    Histo2DPtr _dalitz_3pi;
    Histo1DPtr _h_KpKm,_h_Kppi0;
    Histo2DPtr _dalitz_KpKmpi;
    Histo1DPtr _h_Kppip,_h_K0pip,_h_K0Kp;
    Histo2DPtr _dalitz_K0Kppim;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2017_I1512302);

}
