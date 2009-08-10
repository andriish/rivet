// -*- C++ -*-
#ifndef RIVET_D0_1996_S3324664_HH
#define RIVET_D0_1996_S3324664_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  class D0_1996_S3324664 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    D0_1996_S3324664();

    /// Factory method
    static Analysis* create() {
      return new D0_1996_S3324664();
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

    AIDA::IHistogram1D *_h_deta;
    BinnedHistogram<double> _h_dphi;
    AIDA::IProfile1D *_h_cosdphi_deta;
    //@}

  };


}

#endif

