// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BELLE_2015_I1330289 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BELLE_2015_I1330289);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      _h_spectrum = bookHisto1D(1, 1, 2);
      _nBottom = 0.;
    }


    void findDecayProducts(const Particle & mother, unsigned int & nK0,
			   unsigned int & nKp, unsigned int & nKm,
			   FourMomentum & ptot) {
      for(const Particle & p : mother.children()) {
        int id = p.pdgId();
        if ( id == PID::KPLUS ) {
	  ++nKp;
	  ptot += p.momentum();
	}
        else if (id == PID::KMINUS ) {
	  ++nKm;
	  ptot += p.momentum();
	}
        else if (id == PID::K0S) {
          ++nK0;
	  ptot += p.momentum();
        }
        else if (id == PID::PI0 || id == PID::PIPLUS || id == PID::PIMINUS) {
	  ptot += p.momentum();
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nK0, nKp, nKm, ptot);
        }
        else
	  ptot += p.momentum();
      }
    }
      
    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Loop over bottoms
      for(const Particle& bottom : apply<UnstableParticles>(event, "UFS").particles()) {
	if(bottom.pdgId()!=-521 && bottom.pdgId()!=-511) continue;
	FourMomentum pgamma(0.,0.,0.,0.);
	unsigned int ngamma = 0;
	bool fs = true;
	foreach(const Particle & child, bottom.children()) {
	  if(child.pdgId()==bottom.pdgId()) {
	    fs = false;
	    break;
	  }
	  else if(child.pdgId()==22) {
	    ngamma += 1;
	    pgamma += child.momentum();
	  }
	}
	if(!fs) continue;
	_nBottom += event.weight();
	if(ngamma!=1) continue;
        unsigned int nK0(0),nKp(0),nKm(0);
      	FourMomentum p_tot(0,0,0,0);
        findDecayProducts(bottom, nK0, nKp, nKm, p_tot);
	unsigned int nk = nKp-nKm+nK0;
	if(nk%2==1) {
	  p_tot -= pgamma;
	  _h_spectrum->fill(p_tot.mass()/GeV,event.weight());
	}
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_spectrum, 1e6/_nBottom);
      // multiply by the bin width
      for (unsigned int ix=0;ix<_h_spectrum->numBins();++ix) {
      	_h_spectrum->bins()[ix].scaleW(_h_spectrum->bins()[ix].xWidth());
      }
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_spectrum;
    double _nBottom;
    //@}
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BELLE_2015_I1330289);


}
