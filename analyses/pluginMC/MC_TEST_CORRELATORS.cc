// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Tools/Correlators.hh"
#include "YODA/Scatter2D.h"

namespace Rivet {


  class MC_TEST_CORRELATORS : public CumulantAnalysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    MC_TEST_CORRELATORS() : CumulantAnalysis("MC_TEST_CORRELATORS") {
    }
    //@}

  public:

    /// @name Analysis methods
    //@{
    /// Book histograms and initialise projections before the run
    void init() {
      ChargedFinalState cfs(Cuts::abseta < 1.0);
      declare(cfs, "CFS");
      ChargedFinalState pp(Cuts::abseta < 2.0);
      declare(pp, "PP");
      book(h_c22, "c22", {0.,1.,2.,5.,20.,30.,50.,70.,100});
      ec22 = bookECorrelator<2,2>("ec22",{0.,1.,2.,5.,20.,30.,50.,70.,100});
      pair<int, int> max = getMaxValues(); 
      // Declare correlator projections.
      declare(Correlators(pp, max.first, max.second),"CRS");
    }
    
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const Correlators& c = apply<Correlators>(event,"CRS");
      ec22->fill(apply<ChargedFinalState>(event,"CFS").particles().size(), c);
    }
    
    /// Normalise histograms etc., after the run
    void finalize() {
      stream();
      cnTwoInt(h_c22,ec22);
    }

    //@}
  private:


    /// @name Histograms
    //@{
    Scatter2DPtr h_c22;
    ECorrPtr ec22;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_TEST_CORRELATORS);

}
