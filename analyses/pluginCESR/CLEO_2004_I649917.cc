// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D0 -> KS0 pi0 eta
  class CLEO_2004_I649917 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2004_I649917);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_Kpi  ,1,1,1);
      book(_h_etapi,1,1,2);
      book(_h_Keta ,1,1,3);
      book(_dalitz, "dalitz",50,0.4,1.9,50,0.3,1.8);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & K0 , Particles & pi0,
			   Particles & eta) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::KPLUS || id == PID::KMINUS || id == PID::K0L ||
	     id == PID::PIPLUS || id == PID::PIMINUS) {
	  ++nstable;
	}
	else if (id == PID::K0S) {
	  K0.push_back(p);
	  ++nstable;
	}
	else if (id == PID::PI0) {
	  pi0.push_back(p);
	  ++nstable;
	}
	else if (id == PID::ETA) {
	  eta.push_back(p);
          ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, K0, pi0, eta);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 421 )) {
	unsigned int nstable(0);
	Particles K0, pi0, eta;
	findDecayProducts(meson, nstable, K0, pi0, eta);
	if(nstable !=3) continue;
	if (eta.size()==1&&pi0.size()==1&&K0.size()==1) {
	  double metapi = (eta[0].momentum()+pi0[0].momentum()).mass2();
	  double metaK  = (eta[0].momentum()+ K0[0].momentum()).mass2();
	  double mKpi   = ( K0[0].momentum()+pi0[0].momentum()).mass2();
	  _dalitz ->fill(metapi,mKpi);
	  _h_etapi->fill(metapi);
	  _h_Keta ->fill(metaK );
	  _h_Kpi  ->fill(mKpi);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Kpi  );
      normalize(_h_etapi);
      normalize(_h_Keta );
      normalize(_dalitz );
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kpi,_h_etapi,_h_Keta;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_2004_I649917);

}
