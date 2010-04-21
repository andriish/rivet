// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief MC validation analysis for Z + jets events
  /// @todo More! This analysis just checks the \f$ \eta \f$ distribution at the moment.
  class MC_TTBAR : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    MC_TTBAR()
      : Analysis("MC_TTBAR")
    {
      //setNeedsCrossSection(false);
    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      const ChargedFinalState cfs(-5.0, 5.0);
      addProjection(cfs, "CFS");

      /// @todo Book histograms here, e.g.:
      _hist_nch_eta = bookHistogram1D("nch-eta", 20, -5.0, 5.0);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(event, "CFS");

      foreach (const Particle& p, cfs.particles()) {
        double eta = p.momentum().pseudorapidity();
        _hist_nch_eta->fill(eta, weight);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_hist_nch_eta, 1.0/sumOfWeights());
    }


  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* _hist_nch_eta;
    //@}

  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<MC_TTBAR> plugin_MC_TTBAR;


}
