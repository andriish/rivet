// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief eta' -> pi+pi- gamma decays
  class CRYSTAL_BARREL_1997_I456942 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CRYSTAL_BARREL_1997_I456942);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // Book histograms
      book(_h_m, 1, 1, 1);
    }
    
    void findDecayProducts(const Particle & mother, unsigned int & nstable, unsigned int & ngamma, 
                           unsigned int & npip, unsigned int & npim, FourMomentum & ptot) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if (id == PID::PIMINUS ) {
	  ++npim;
          ++nstable;
	  ptot += p.momentum();
	}
        else if (id == PID::PIPLUS) {
          ++npip;
          ++nstable;
	  ptot += p.momentum();
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable, ngamma,npip,npim,ptot);
        }
        else if (id == PID::GAMMA) {
	  ++ngamma;
          ++nstable;
        }
        else
          ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Loop over eta' mesons
      for (const Particle& p :  apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==331)) {
	unsigned nstable(0),ngamma(0),npip(0),npim(0);
	FourMomentum ptot;
	findDecayProducts(p,nstable,ngamma,npip,npim,ptot);
	if(nstable==3 && npim==1 && npip==1 && ngamma==1)
	  _h_m->fill(ptot.mass()/MeV);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_m,1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_m;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CRYSTAL_BARREL_1997_I456942);

}
