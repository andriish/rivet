// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Jet masses for 14, 22, 34.8 and 43,5 GeV
  class TASSO_1989_I279165 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TASSO_1989_I279165);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      const Thrust thrust(cfs);
      declare(thrust, "Thrust");
      declare(Hemispheres(thrust), "Hemispheres");

      // Compute histogram indexes from beam energy
      int offset = 0;
      if      (beamEnergyMatch(14.0*GeV)) offset=1;
      else if (beamEnergyMatch(22.0*GeV)) offset=2;
      else if (beamEnergyMatch(34.8*GeV)) offset=3;
      else if (beamEnergyMatch(43.5*GeV)) offset=4;
      else MSG_ERROR("Beam energy " << sqrtS() << " not supported!");

      // Book histograms
      book(_h_diff , 1, 1, offset);
      book(_h_heavy, 2, 1, offset);
      book(_h_light, 3, 1, offset);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");

      if (cfs.particles().size() < 3 ) vetoEvent;

      const Hemispheres& hemi = apply<Hemispheres>(event, "Hemispheres");
      _h_heavy->fill(hemi.scaledM2high());
      _h_light->fill(hemi.scaledM2low() );
      _h_diff ->fill(hemi.scaledM2diff());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_diff );
      normalize(_h_heavy);
      normalize(_h_light);
    }
    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_diff,_h_heavy,_h_light;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TASSO_1989_I279165);


}
