// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BESIII_2019_I1724880 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BESIII_2019_I1724880);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      _h_n = bookHisto1D(1, 1, 1);

    }
    
    void findChildren(const Particle & p,int & nCharged) {
      foreach(const Particle &child, p.children()) {
	if(child.children().empty()) {
	  if(PID::isCharged(child.pdgId())) ++nCharged;
	}
	else
	  findChildren(child,nCharged);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      foreach (const Particle& p, apply<FinalState>(event, "UFS").particles(Cuts::pid==441)) {
	int nCharged(0);
	findChildren(p,nCharged);
	_h_n->fill(min(nCharged,8),event.weight());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_n,2.);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_n;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BESIII_2019_I1724880);


}