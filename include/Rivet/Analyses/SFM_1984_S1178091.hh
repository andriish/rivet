// -*- C++ -*-
#ifndef RIVET_SFM_1984_S1178091_HH
#define RIVET_SFM_1984_S1178091_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  class SFM_1984_S1178091 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    SFM_1984_S1178091();

    /// Factory method
    static Analysis* create() {
      return new SFM_1984_S1178091();
    }
    //@}


  public:

    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

  private:

    // Data members like post-cuts event weight counters go here

  private:

    /// @name Histograms
    //@{

    AIDA::IHistogram1D *_hist_multiplicity_inel_30;
    AIDA::IHistogram1D *_hist_multiplicity_inel_45;
    AIDA::IHistogram1D *_hist_multiplicity_inel_53;
    AIDA::IHistogram1D *_hist_multiplicity_inel_63;
    AIDA::IHistogram1D *_hist_multiplicity_nsd_30;
    AIDA::IHistogram1D *_hist_multiplicity_nsd_45;
    AIDA::IHistogram1D *_hist_multiplicity_nsd_53;
    AIDA::IHistogram1D *_hist_multiplicity_nsd_63;
    //@}

  };


}

#endif

