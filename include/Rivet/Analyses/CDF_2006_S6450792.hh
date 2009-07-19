// -*- C++ -*-
#ifndef RIVET_CDF_2006_S6450792_HH
#define RIVET_CDF_2006_S6450792_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  class CDF_2006_S6450792 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CDF_2006_S6450792();

    /// Factory method
    static Analysis* create() {
      return new CDF_2006_S6450792();
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

    /// @name Histograms
    //@{

    AIDA::IHistogram1D *_h_jet_pt;
    //@}

  };


}

#endif

