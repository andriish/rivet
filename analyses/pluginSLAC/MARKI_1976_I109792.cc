// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Charged particle spectra between 3.0 and 7.4 GeV
  class MARKI_1976_I109792 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MARKI_1976_I109792);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState fs;
      declare(fs, "FS");
      unsigned int iloc(0);
      if      (beamEnergyMatch(3.0*GeV)) iloc = 8;
      else if (beamEnergyMatch(4.8*GeV)) iloc = 7;
      else if (beamEnergyMatch(5.8*GeV)) iloc = 6;
      else if (beamEnergyMatch(6.2*GeV)) iloc = 5;
      else if (beamEnergyMatch(6.6*GeV)) iloc = 4;
      else if (beamEnergyMatch(7.0*GeV)) iloc = 3;
      else if (beamEnergyMatch(7.4*GeV)) iloc = 2;
      assert(iloc!=0);
      book(_h_x, iloc   ,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& fs = apply<ChargedFinalState>(event, "FS");
      if(fs.particles().size()==2 &&
         fs.particles()[0].abspid()==PID::MUON &&
         fs.particles()[1].abspid()==PID::MUON) vetoEvent;
      for (const Particle& p : fs.particles()) {
        const Vector3 mom3 = p.p3();
        double pp = mom3.mod();
        double x = 2.*pp/sqrtS();
        _h_x->fill(x);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_x,crossSection()*sqr(sqrtS())/sumOfWeights()/microbarn);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_x;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MARKI_1976_I109792);


}
