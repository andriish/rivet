// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D+ -> K- pi+ pi+
  class CLEOC_2008_I780363 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEOC_2008_I780363);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      
      book(_h_Kpiall,1,1,1);
      book(_h_pipi  ,1,2,1);

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
	else if (id == PID::K0S||id == PID::K0L) {
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
      // parameters from the efficiency function, table 1 in paper
      static const double E1   = -0.0153;
      static const double E2   = -0.030;
      static const double E3   =  0.162;
      static const double Exy  = -0.053;
      static const double Exyn =  0.673;
      // static const double Eth[3] = {4.25,4.25,2.907};
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 411)) {
	unsigned int nstable(0);
	Particles pip, pim, pi0, Kp , Km, K0;
	findDecayProducts(meson, nstable, pip, pim, pi0, Kp , Km, K0);
	if(meson.pid()<0) {
	  swap(pim,pip);
	  swap(Kp,Km);
	}
	if(pip.size()==2&&Km.size()==1) {
	  // kinematic variables
	  double x[3]  = {(Km[0].momentum() +pip[0].momentum()).mass2(),
			  (Km[0].momentum() +pip[1].momentum()).mass2(),
			  (pip[0].momentum()+pip[1].momentum()).mass2()};
	  if(x[1]<x[0]) swap(x[0],x[1]);
	  // calculate the efficiency, eqns 6,7 from paper
	  // double xmax[3] = {sqr(meson.mass()-pip[0].mass()),
	  // 		    sqr(meson.mass()-pip[0].mass()),
	  // 		    sqr(meson.mass()-Km[0].mass())};
	  double xh = x[0]-1.5,yh=x[1]-1.5;
	  double eff = (1.+E1*(xh+yh)+E2*(sqr(xh)+sqr(yh))+E3*(pow(xh,3)+pow(yh,3))
			+Exy*xh*yh+Exyn*xh*yh*(xh+yh));
	  // double T=1.;
	  // for(unsigned int ix=2;ix<3;++ix) {
	  //   double arg = Eth[ix]*abs(x[ix]-xmax[ix]);
	  //   if(arg<0.5*M_PI) T *=sin(arg);
	  // }
	  // eff *=T;
	  // fill plots
	  _h_Kpiall->fill( x[0],eff);
	  _h_Kpiall->fill( x[1],eff);
	  _h_pipi  ->fill( x[2],eff);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Kpiall);
      normalize(_h_pipi  );
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_Kpiall, _h_pipi;
    //@}


  };


  DECLARE_RIVET_PLUGIN(CLEOC_2008_I780363);

}
