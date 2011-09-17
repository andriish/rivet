// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "HepMC/HepMCDefs.h"

namespace Rivet {


  /// @brief Analysis for the generated cross section
  class MC_XS : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    MC_XS()
      : Analysis("MC_XS")
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      _h_XS = bookDataPointSet("XS", 1, 0.0, 1.0);
      _h_N = bookHistogram1D("N", 1, 0.0, 1.0);
      _h_pmXS = bookHistogram1D("pmXS", 2, -1.0, 1.0);
      _h_pmN = bookHistogram1D("pmN", 2, -1.0, 1.0);
      _mc_xs=_mc_error=0.;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      _h_N->fill(0.5,1.);
      _h_pmXS->fill(0.5*(event.weight()>0?1.:-1),abs(event.weight()));
      _h_pmN->fill(0.5*(event.weight()>0?1.:-1),1.);
#ifdef HEPMC_HAS_CROSS_SECTION
      _mc_xs    = event.genEvent().cross_section()->cross_section();
      _mc_error = event.genEvent().cross_section()->cross_section_error();
#endif
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_pmXS, crossSection()/sumOfWeights());
#ifndef HEPMC_HAS_CROSS_SECTION
      _mc_xs=crossSection();
      _mc_error=0.0;
#endif
      std::vector<double> xs,err;
      xs.push_back(_mc_xs);
      err.push_back(_mc_error);
      _h_XS->setCoordinate(1,xs,err);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    AIDA::IDataPointSet * _h_XS;
    AIDA::IHistogram1D *  _h_N;
    AIDA::IHistogram1D *  _h_pmXS;
    AIDA::IHistogram1D *  _h_pmN;
    double _mc_xs, _mc_error;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_XS);

}
