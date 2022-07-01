// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D0 -> K- pi+ eta
  class BELLE_2020_I1785816 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2020_I1785816);


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
      book(_dalitz, "dalitz",50,0.3,1.8,50,0.4,1.9);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip, Particles & pim ,
			   Particles & Kp, Particles & Km, Particles & eta) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if (id == PID::PI0|| id == PID::K0S || id == PID::K0L) {
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
	else if (id == PID::ETA) {
	  eta.push_back(p);
          ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, Kp, Km, eta);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 421 )) {
	unsigned int nstable(0);
	Particles pip, pim, Kp, Km, eta;
	findDecayProducts(meson, nstable, pip, pim, Kp, Km, eta);
	if(nstable !=3) continue;
	if(meson.pid()<0) {
	  swap(pim,pip);
	  swap(Km ,Kp );
	}
	if (eta.size()==1&&pip.size()==1&&Km.size()==1) {
	  double metapi = (eta[0].momentum()+pip[0].momentum()).mass2();
	  double metaK  = (eta[0].momentum()+ Km[0].momentum()).mass2();
	  double mKpi   = ( Km[0].momentum()+pip[0].momentum()).mass2();
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


  RIVET_DECLARE_PLUGIN(BELLE_2020_I1785816);

}
