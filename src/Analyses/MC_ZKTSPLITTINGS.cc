// -*- C++ -*-
#include "Rivet/Analyses/MC_JetSplittings.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetYODA.hh"

namespace Rivet {


  /// @brief MC validation analysis for Z + jets events
  class MC_ZKTSPLITTINGS : public MC_JetSplittings {
  public:

    /// Default constructor
    MC_ZKTSPLITTINGS()
      : MC_JetSplittings("MC_ZKTSPLITTINGS", 4, "Jets")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      FinalState fs;
      ZFinder zfinder(fs, -3.5, 3.5, 25.0*GeV, ELECTRON, 65.0*GeV, 115.0*GeV, 0.2, true, true);
      addProjection(zfinder, "ZFinder");
      FastJets jetpro(zfinder.remainingFinalState(), FastJets::KT, 0.6);
      addProjection(jetpro, "Jets");

      MC_JetSplittings::init();
    }



    /// Do the analysis
    void analyze(const Event & e) {
      const ZFinder& zfinder = applyProjection<ZFinder>(e, "ZFinder");
      if (zfinder.bosons().size()!=1) {
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
  DECLARE_RIVET_PLUGIN(MC_ZKTSPLITTINGS);

}
