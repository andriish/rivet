// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief eta' -> gamma gamma pi0
  class BESIII_2016_I1504943 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2016_I1504943);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      for(unsigned int ix=0;ix<3;++ix)
	book(_h_br[ix],1,1,1+ix);
      book(_h_m, 2, 1, 1);
      book(_netap, "TMP/netap");
    }
    
    void findDecayProducts(const Particle & mother, unsigned int & nstable, unsigned int & ngamma, 
                           unsigned int & npi0, FourMomentum & pgamma,
			   bool &omega, bool &nr) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if (id == PID::PI0 ) {
	  ++npi0;
          ++nstable;
	}
        else if (id == PID::GAMMA) {
          ++ngamma;
          ++nstable;
	  pgamma += p.momentum();
        }
        else if ( !p.children().empty() ) {
	  if (p.pid()==223) omega=true;
	  nr = false;
	  findDecayProducts(p, nstable, ngamma,npi0,pgamma,omega,nr);
        }
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Loop over eta' mesons
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==331)) {
	unsigned nstable(0),ngamma(0),npi0(0);
	bool omega(false),nr(true);
	FourMomentum pgamma;
	findDecayProducts(p,nstable,ngamma,npi0,pgamma,omega,nr);
	_netap->fill();
	if(nstable==3 && npi0==1 && ngamma==2) {
	  _h_m->fill(pgamma.mass2());
	  _h_br[0]->fill(.5);
	  if(omega) _h_br[1]->fill(.5);
	  if(nr) _h_br[2]->fill(0.5);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // eta' width in kev
      double gammaEtap = 0.188e3;
      scale(_h_m, gammaEtap/ *_netap);
      for(unsigned int ix=0;ix<3;++ix)
	scale(_h_br[ix], 1e4/ *_netap);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_m,_h_br[3];
    CounterPtr _netap;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2016_I1504943);

}
