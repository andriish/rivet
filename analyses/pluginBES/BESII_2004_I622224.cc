// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BESII_2004_I622224 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BESII_2004_I622224);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState fs;
      declare(fs, "FS");
      unsigned int iloc(0);
      if(fuzzyEquals(sqrtS()/GeV, 2.2 , 1E-3))
	iloc = 1;
      else if(fuzzyEquals(sqrtS()/GeV, 2.6 , 1E-3))
	iloc = 2;
      else if(fuzzyEquals(sqrtS()/GeV, 3.0 , 1E-3))
	iloc = 3;
      else if(fuzzyEquals(sqrtS()/GeV, 3.2 , 1E-3))
	iloc = 4;
      else if(fuzzyEquals(sqrtS()/GeV, 4.6 , 1E-3))
	iloc = 5;
      else if(fuzzyEquals(sqrtS()/GeV, 4.8 , 1E-3))
	iloc = 6;
      assert(iloc!=0);
      _h_ln =  bookHisto1D( iloc   ,1,1);
      _h_weight = 0.;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& fs = apply<ChargedFinalState>(event, "FS");
      if(fs.particles().size()==2 &&
	 abs(fs.particles()[0].pdgId())==13 &&
	 abs(fs.particles()[1].pdgId())==13) vetoEvent;
      foreach (const Particle& p, fs.particles()) {
	const Vector3 mom3 = p.p3();
	double pp = mom3.mod();
	double xi = -log(2.*pp/sqrtS());
	_h_ln->fill(xi,event.weight());
      }
      _h_weight+= event.weight();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_ln,1./_h_weight);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_ln;
    double _h_weight;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BESII_2004_I622224);


}
