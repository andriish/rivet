// -*- C++ -*-
#include "Rivet/Analysis.hh"

#ifndef ENABLE_HEPMC_3
#include "HepMC/HepMCDefs.h"
#endif

namespace Rivet {

  /// @brief Analysis for the generated cross section
  class MC_XS : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_XS);

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      /// @todo Convert to Scatter1D or Counter
      book(_h_XS, "XS");
      book(_h_N, "N");
      //book(_h_pmXS, "pmXS", 2, -1.0, 1.0);
      //book(_h_pmN, "pmN", 2, -1.0, 1.0);
      book(_c_N, "_aux_N");

      _mc_xs = _mc_error = 0.;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      _c_N->fill();
      //_h_pmXS->fill(0.5);
      //_h_pmN ->fill(0.5);
      
      #if defined ENABLE_HEPMC_3
      //@todo HepMC3::GenCrossSection methods aren't const accessible :(
      RivetHepMC::GenCrossSection gcs = *(event.genEvent()->cross_section());
      _mc_xs    = gcs.xsec();
      _mc_error = gcs.xsec_err();
      #elif defined HEPMC_HAS_CROSS_SECTION
      _mc_xs    = event.genEvent()->cross_section()->cross_section();
      _mc_error = event.genEvent()->cross_section()->cross_section_error();
      #endif // VERSION_CODE >= 3000000
      
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      //scale(_h_pmXS, crossSection()/sumOfWeights());
      const double N = _c_N->numEntries();
      _h_N->addPoint(0.5, N, 0.5, sqrt(N));
      #ifndef HEPMC_HAS_CROSS_SECTION
      _mc_xs = crossSection();
      _mc_error = 0.0;
      #endif
      _h_XS->addPoint(0, _mc_xs, 0.5, _mc_error);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    CounterPtr _c_N;
    Scatter2DPtr _h_XS;
    Scatter2DPtr _h_N;
    //Histo1DPtr _h_pmXS;
    //Histo1DPtr _h_pmN;
    double _mc_xs, _mc_error;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_XS);

}
