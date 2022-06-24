// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D+ -> K+K-pi+
  class CLEO_2008_I791716 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2008_I791716);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      book(_h_Kppi,1,1,1);
      book(_h_Kmpi,1,1,2);
      book(_h_KK  ,1,1,3);
      book(_dalitz, "dalitz",50,0.,2.,50,0.9,3.1);
    }


    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim ,
			   Particles & Kp  , Particles & Km ) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::KPLUS) {
	  Kp.push_back(p);
	  ++nstable;
	}
	else if( id == PID::KMINUS) {
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
	else if (id == PID::K0S || id == PID::PI0 || id == PID::K0L) {
	  ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, Kp, Km);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").
	    particles(Cuts::abspid== 411)) {
	unsigned int nstable(0);
	Particles pip, pim, Kp, Km;
	findDecayProducts(meson, nstable, pip, pim, Kp, Km);
	if(nstable !=3) continue;
	if(meson.pid()<0) {
	  swap(pim,pip);
	  swap(Km ,Kp );
	}
	if(Km.size()==1&&Kp.size()==1&&pip.size()==1) {
	  double mminus = (Km[0].momentum()+pip[0].momentum() ).mass2();
	  double mplus  = (Kp[0].momentum()+pip[0].momentum() ).mass2();
	  double mKK    = (Kp[0].momentum()+Km[0].momentum()).mass2();
	  _h_Kppi->fill(mplus);
	  _h_Kmpi->fill(mminus);
	  _h_KK  ->fill(mKK);
	  _dalitz->fill(mminus,mKK); 
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Kppi);
      normalize(_h_Kmpi);
      normalize(_h_KK  );
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kppi,_h_Kmpi,_h_KK;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_2008_I791716);

}
