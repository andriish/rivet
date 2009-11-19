// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
/// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...

namespace Rivet {


  class MC_LHC_TTBAR : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    MC_LHC_TTBAR()
      : Analysis("MC_LHC_TTBAR")
    {
      /// @todo Set approriate for your analysis
      setBeams(PROTON, PROTON);
   
      /// @todo Set whether your finalize method needs the generator cross section
      setNeedsCrossSection(false);

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
       _histPseudorapidity = bookHistogram1D("pseudorap", 20, -5.0, 5.0);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(event, "CFS");
   
      foreach (const Particle& p, cfs.particles()) {
        double eta = p.momentum().pseudorapidity();
        _histPseudorapidity->fill(eta, weight);
     }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
    scale(_histPseudorapidity, 1.0/sumOfWeights());
   
    }


  private:

    /// @name Histograms
    //@{

    AIDA::IHistogram1D *_histPseudorapidity;
    //@}

  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<MC_LHC_TTBAR> plugin_MC_LHC_TTBAR;


}
