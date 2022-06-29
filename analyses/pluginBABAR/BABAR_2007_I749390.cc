// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D0 -> K+K-pi0
  class BABAR_2007_I749390 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2007_I749390);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_Kppi ,1,1,1);
      book(_h_Kmpi ,1,1,2);
      book(_h_KK   ,1,1,3);
      book(_dalitz, "dalitz",50,0.3,2.1,50,0.3,2.1);
    }


    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pi0 , Particles & Kp, Particles & Km) {
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
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 421 )) {
	unsigned int nstable(0);
	Particles pi0, Kp , Km;
	findDecayProducts(meson, nstable, pi0, Kp , Km);
	if(nstable !=3) continue;
	if(meson.pid()<0) {
	  swap(Kp,Km);
	}
	if (Km.size()==1&&Kp.size()==1&&pi0.size()==1) {
	  double mplus  = (Kp[0].momentum()+pi0[0].momentum()).mass2();
	  double mminus = (Km[0].momentum()+pi0[0].momentum()).mass2();
	  double mKK    = (Kp[0].momentum()+Km [0].momentum()).mass2();
	  _h_KK  ->fill(mKK);
	  _h_Kppi->fill(mplus);
	  _h_Kmpi->fill(mminus);
	  _dalitz->fill(mminus,mplus);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_KK  );
      normalize(_h_Kmpi);
      normalize(_h_Kppi);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_KK,_h_Kmpi,_h_Kppi;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2007_I749390);

}
