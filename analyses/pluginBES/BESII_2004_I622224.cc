// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Charged particle spectra between 2.2 and 4.8 GeV
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
      unsigned int iloc = 0;
      if      (beamEnergyMatch(2.2*GeV)) iloc = 1;
      else if (beamEnergyMatch(2.6*GeV)) iloc = 2;
      else if (beamEnergyMatch(3.0*GeV)) iloc = 3;
      else if (beamEnergyMatch(3.2*GeV)) iloc = 4;
      else if (beamEnergyMatch(4.6*GeV)) iloc = 5;
      else if (beamEnergyMatch(4.8*GeV)) iloc = 6;
      else MSG_ERROR("Beam energy not supported!");
      assert(iloc != 0);
      book(_h_ln, iloc, 1, 1);
      book(_h_weight, "TMP/Weight");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& fs = apply<ChargedFinalState>(event, "FS");
      if(fs.particles().size()==2 &&
         abs(fs.particles()[0].pid())==13 &&
         abs(fs.particles()[1].pid())==13) vetoEvent;
      for (const Particle& p : fs.particles()) {
        const Vector3 mom3 = p.p3();
        double pp = mom3.mod();
        double xi = -log(2.*pp/sqrtS());
        _h_ln->fill(xi);
      }
      _h_weight->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_ln,1./_h_weight->sumW());
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_ln;
    CounterPtr _h_weight;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BESII_2004_I622224);


}
