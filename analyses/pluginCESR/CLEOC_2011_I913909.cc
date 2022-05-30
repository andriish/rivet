// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D0 -> pi+pi-pi0
  class CLEOC_2011_I913909 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOC_2011_I913909);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_h_Kpi ,1,1,1);
      book(_h_pipi,1,1,2);
      book(_dalitz, "dalitz",50,0.,1.9,50,0.3,3.1);
    }


    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & K0, Particles & pi0) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::KPLUS  || id == PID::KMINUS ||
	     id == PID::PIPLUS || id == PID::PIMINUS || id == PID::K0L) {
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
	  findDecayProducts(p, nstable, K0, pi0);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid== 421 )) {
	unsigned int nstable(0);
	Particles K0, pi0;
	findDecayProducts(meson, nstable, K0, pi0);
	if(nstable !=3) continue;
	if (K0.size()==1&&pi0.size()==2) {
	  double mpipi = (pi0[0].momentum()+pi0[1].momentum()).mass2();
	  _h_pipi->fill(mpipi);
	  for(unsigned int ix=0;ix<2;++ix) {
	    double mKpi  = (K0[0].momentum()+pi0[ix].momentum()).mass2();
	    _h_Kpi->fill(mKpi);
	    _dalitz->fill(mpipi,mKpi);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pipi);
      normalize(_h_Kpi);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kpi, _h_pipi;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOC_2011_I913909);

}
