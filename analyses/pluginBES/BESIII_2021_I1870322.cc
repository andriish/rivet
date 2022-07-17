// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  D_s+ -> pi+ pi+ pi- eta
  class BESIII_2021_I1870322 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1870322);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_pippip   ,1,1,1);
      book(_h_pippim   ,1,1,2);
      book(_h_pipeta   ,1,1,3);
      book(_h_pimeta   ,1,1,4);
      book(_h_3pi      ,1,1,5);
      book(_h_pippimeta,1,1,6);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , Particles & eta) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::KPLUS || id == PID::KMINUS ||
	     id == PID::K0S||id == PID::K0L || id == PID::PI0) {
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
	else if (id == PID::ETA) {
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
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 431 )) {
	unsigned int nstable(0);
	Particles pip, pim, eta;
	findDecayProducts(meson, nstable, pip, pim, eta);
	if(nstable !=4) continue;
	if(meson.pid()<0) {
	  swap(pim,pip);
	}
	if (eta.size()==1&&pim.size()==1&&pip.size()==2) {
	  double metapipi[2] = {(eta[0].momentum()+pim[0].momentum()+pip[0].momentum()).mass(),
				(eta[0].momentum()+pim[0].momentum()+pip[1].momentum()).mass()};
	  if(metapipi[0]<GeV || metapipi[1]<GeV) continue;
	  _h_pippip->fill((pip[0].momentum()+pip[1].momentum()).mass());
	  _h_pippim->fill((pim[0].momentum()+pip[0].momentum()).mass());
	  _h_pippim->fill((pim[0].momentum()+pip[1].momentum()).mass());
	  _h_pipeta->fill((eta[0].momentum()+pip[0].momentum()).mass());
	  _h_pipeta->fill((eta[0].momentum()+pip[1].momentum()).mass());
	  _h_pimeta->fill((eta[0].momentum()+pim[0].momentum()).mass());
	  _h_3pi   ->fill((pim[0].momentum()+pip[0].momentum()+pip[1].momentum()).mass());
	  _h_pippimeta->fill(metapipi[0]);
	  _h_pippimeta->fill(metapipi[1]);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pippip   );
      normalize(_h_pippim   );
      normalize(_h_pipeta   );
      normalize(_h_pimeta   );
      normalize(_h_3pi      );
      normalize(_h_pippimeta);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pippip,_h_pippim,_h_pipeta,_h_pimeta,_h_3pi,_h_pippimeta;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1870322);

}
