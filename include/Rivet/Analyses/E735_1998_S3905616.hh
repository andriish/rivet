// -*- C++ -*-
#ifndef RIVET_E735_1998_S3905616_HH
#define RIVET_E735_1998_S3905616_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  class E735_1998_S3905616 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    E735_1998_S3905616();

    /// Factory method
    static Analysis* create() {
      return new E735_1998_S3905616();
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

    AIDA::IHistogram1D *_hist_multiplicity;
    //@}

  };


}

#endif

