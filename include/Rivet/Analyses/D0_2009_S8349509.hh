// -*- C++ -*-
#ifndef RIVET_D0_2009_S8349509_HH
#define RIVET_D0_2009_S8349509_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  class D0_2009_S8349509 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    D0_2009_S8349509();

    /// Factory method
    static Analysis* create() {
      return new D0_2009_S8349509();
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
    AIDA::IHistogram1D *_h_dphi_jet_Z25;
    AIDA::IHistogram1D *_h_dphi_jet_Z45;

    AIDA::IHistogram1D *_h_dy_jet_Z25;
    AIDA::IHistogram1D *_h_dy_jet_Z45;

    AIDA::IHistogram1D *_h_yboost_jet_Z25;
    AIDA::IHistogram1D *_h_yboost_jet_Z45;
    //@}
    
    double _inclusive_Z_sumofweights;

  };


}

#endif

