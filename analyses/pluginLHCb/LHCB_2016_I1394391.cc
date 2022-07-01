// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  D0 -> KS) K+/- pi-/+
  class LHCB_2016_I1394391 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2016_I1394391);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_Kmpip,1,1,1);
      book(_h_K0pip,1,1,2);
      book(_h_K0Km ,1,1,3);
      book(_h_Kppim,2,1,1);
      book(_h_K0pim,2,1,2);
      book(_h_K0Kp ,2,1,3);
      book(_dalitz [0],"dalitz_1",50,0.3,2.0,50,0.3,2.);
      book(_dalitz [1],"dalitz_2",50,0.3,2.0,50,0.3,2.);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , Particles & K0, Particles & Kp,  Particles & Km) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if (id == PID::PI0 || id == PID::K0L) {
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
	else if (id == PID::KPLUS) {
	  Kp.push_back(p);
	  ++nstable;
	}
	else if (id == PID::KMINUS) {
	  Km.push_back(p);
	  ++nstable;
	}
	else if (id == PID::K0S) {
	  K0.push_back(p);
          ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, K0, Kp, Km);
	}
	else
	  ++nstable;
      }
    }

    double efficiency(const double & x, const double &y) {
      double X=x-2., Y=y-1.;
      static const double E0 = 5.8096, Ex = -3.645, Ey = -3.174, Ex2=  0.831,
	Exy = 2.131, Ey2 = 4.43, Ex3 = -0.427, Ex2y = 2.65, Exy2 = 1.50, Ey3 = -3.92;
      return E0 + Ex*X + Ey*Y + Ex2*sqr(X) + Ey2*sqr(Y) + Exy*X*Y +
	Ex3*pow(X,3) + Ex2y*sqr(X)*Y + Exy2*X*sqr(Y) + Ey3*pow(Y,3);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").
	    particles(Cuts::abspid== 421)) {
	unsigned int nstable(0);
	Particles pip, pim, K0,Kp,Km;
	findDecayProducts(meson, nstable, pip, pim, K0,Kp,Km);
	if(nstable !=3 || K0.size()!=1) continue;
	if(meson.pid()<0) {
	  swap(pim,pip);
	  swap(Km,Kp);
	}
	// K0S K- pi+
	if( Km.size()==1 && pip.size()==1) {
	  double mK0pip = (K0[0].momentum()+pip[0].momentum() ).mass2();
	  double mKmpip = (Km[0].momentum()+pip[0].momentum() ).mass2();
	  double mKK    = (K0[0].momentum()+Km [0].momentum() ).mass2();
	  double eff = efficiency(mKK,mK0pip);
	  _h_K0Km ->fill(mKK   ,eff);
	  _h_K0pip->fill(mK0pip,eff);
	  _h_Kmpip->fill(mKmpip,eff);
	  _dalitz[0]->fill(mKmpip,mK0pip); 
	}
	// K0S K+ pi-
	else if( Kp.size()==1 && pim.size()==1) {
	  double mK0pim = (K0[0].momentum()+pim[0].momentum() ).mass2();
	  double mKppim = (Kp[0].momentum()+pim[0].momentum() ).mass2();
	  double mKK    = (K0[0].momentum()+Kp [0].momentum() ).mass2();
	  double eff = efficiency(mKK,mK0pim);
	  _h_K0Kp ->fill(mKK   ,eff);
	  _h_K0pim->fill(mK0pim,eff);
	  _h_Kppim->fill(mKppim,eff);
	  _dalitz[1]->fill(mKppim,mK0pim); 
	}
      }
    }


    /// Normalise histograms etc., after the runbook
    void finalize() {
      normalize(_h_Kmpip);
      normalize(_h_K0pip);
      normalize(_h_K0Km );
      normalize(_h_Kppim);
      normalize(_h_K0pim);
      normalize(_h_K0Kp );
      normalize(_dalitz [0]);
      normalize(_dalitz [1]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kmpip, _h_K0pip, _h_K0Km;
    Histo1DPtr _h_Kppim, _h_K0pim, _h_K0Kp;
    Histo2DPtr _dalitz[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2016_I1394391);

}
