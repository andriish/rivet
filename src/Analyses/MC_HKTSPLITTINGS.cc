// -*- C++ -*-
#include "Rivet/Analyses/MC_JetSplittings.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetYODA.hh"

namespace Rivet {


  /// @brief MC validation analysis for higgs [-> tau tau] + jets events
  class MC_HKTSPLITTINGS : public MC_JetSplittings {
  public:

    /// Default constructor
    MC_HKTSPLITTINGS()
      : MC_JetSplittings("MC_HKTSPLITTINGS", 4, "Jets")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      FinalState fs;
      ZFinder hfinder(fs, -3.5, 3.5, 25.0*GeV, TAU, 115.0*GeV, 125.0*GeV, 0.0, false, false);
      addProjection(hfinder, "Hfinder");
      FastJets jetpro(hfinder.remainingFinalState(), FastJets::KT, 0.6);
      addProjection(jetpro, "Jets");

      MC_JetSplittings::init();
    }



    /// Do the analysis
    void analyze(const Event & e) {
      const ZFinder& hfinder = applyProjection<ZFinder>(e, "Hfinder");
      if (hfinder.bosons().size()!=1) {
        vetoEvent;
      }

      MC_JetSplittings::analyze(e);
    }


    /// Finalize
    void finalize() {
      MC_JetSplittings::finalize();
    }

    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_HKTSPLITTINGS);

}
