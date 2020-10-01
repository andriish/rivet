// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// Hadronic charged-multiplicity measurement between 14 and 43.6 GeV
  class TASSO_1989_I277658 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TASSO_1989_I277658);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Projections
      const ChargedFinalState cfs;
      declare(cfs, "CFS");

      // Compute histo index offsets from beam energies
      int offset = 0;
      if      (beamEnergyMatch(14.0*GeV)) offset = 1;
      else if (beamEnergyMatch(22.0*GeV)) offset = 2;
      else if (beamEnergyMatch(34.8*GeV)) offset = 3;
      else if (beamEnergyMatch(43.6*GeV)) offset = 4;
      else MSG_WARNING("CoM energy of events sqrt(s) = " << sqrtS()/GeV
                       << " doesn't match any available analysis energy .");

      // Book histograms
      book(_histCh, 5, 1, offset);
      book(_histTotal, 2, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      MSG_DEBUG("Total charged multiplicity = " << cfs.size());
      _histCh->fill(cfs.size());
      _histTotal->fill(sqrtS(),cfs.size());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_histCh, 2.0/sumOfWeights()); // bin width (2)
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _histCh;
    Profile1DPtr _histTotal;
    //@}
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TASSO_1989_I277658);


}
