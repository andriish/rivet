// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Tools/Correlators.hh"


namespace Rivet {


  class TEST : public CumulantAnalysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    TEST() : CumulantAnalysis("TEST") {
    }
    //@}

  public:

    /// @name Analysis methods
    //@{
    /// Book histograms and initialise projections before the run
    void init() {

      ChargedFinalState cfs(-1.0, 1.0);
      declare(cfs, "CFS");
      ChargedFinalState pp(Cuts::abseta < 2.0);
      declare(pp, "PP");
      h_c22 = bookScatter2D("c22",120,0,120);
      h_c23 = bookScatter2D("c23",120,0,120);
      ec22 = bookECorrelator<2,2>("ec22",h_c22);
      ec23 = bookECorrelator<3,2>("ec32",h_c22);
      pair<int, int> max = getMaxValues(); 
      // Declare correlator projections.
      declare(Correlators(pp, max.first, max.second),"CRS");
    }
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      ec22->fill(apply<ChargedFinalState>(event,"CFS").particles().size(), 
        apply<Correlators>(event,"CRS"), event.weight());
      ec23->fill(apply<ChargedFinalState>(event,"CFS").particles().size(), 
        apply<Correlators>(event,"CRS"), event.weight());
    }
    /// Normalise histograms etc., after the run
    void finalize() {
      CumulantAnalysis::finalize();
      cnTwoInt(h_c22,ec22);
    }


    //@}
  private:


    /// @name Histograms
    //@{
    Scatter2DPtr h_c22;
    ECorrPtr ec22;
    Scatter2DPtr h_c23;
    ECorrPtr ec23;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TEST);

}
