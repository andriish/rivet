// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D0 -> and K+K-
  class BESIII_2020_I1799437 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2020_I1799437);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // Histograms
      book(_h_K0Km,1,1,3);
      book(_h_K0Kp,1,1,2);
      book(_h_KpKm,1,1,1);
      book(_dalitz, "dalitz",50,0.9,1.9,50,0.9,1.9);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & K0, Particles & Kp , Particles & Km ) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::PI0    || id == PID::K0L ||
	     id == PID::PIPLUS || id == PID::PIMINUS) {
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
	  findDecayProducts(p, nstable, K0, Kp, Km);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").
	    particles(Cuts::abspid== 421)) {
	unsigned int nstable(0);
	Particles K0,Kp,Km;
	findDecayProducts(meson, nstable, K0, Kp, Km);
	if(nstable !=3) continue;
	if(meson.pid()<0) {
	  swap(Km,Kp);
	}
	if(Km.size()==1&&Kp.size()==1&&K0.size()==1) {
	  double mminus = (Km[0].momentum()+K0[0].momentum() ).mass2();
	  double mplus  = (Kp[0].momentum()+K0[0].momentum() ).mass2();
	  double mKK  = (Kp[0].momentum()+Km[0].momentum()).mass2();
	  _h_K0Kp->fill(mplus);
	  _h_K0Km->fill(mminus);
	  _h_KpKm->fill(mKK);
	  _dalitz->fill(mKK,mplus); 
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_K0Km);
      normalize(_h_K0Kp);
      normalize(_h_KpKm);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_K0Km,_h_KpKm,_h_K0Kp;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2020_I1799437);

}
