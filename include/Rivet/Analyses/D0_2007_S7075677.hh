// -*- C++ -*-
#ifndef RIVET_D0_2007_S7075677_HH
#define RIVET_D0_2007_S7075677_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief Measurement of D0 Run II Z pT diff cross-section shape
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  /// @author Frank Siegert
  class D0_2007_S7075677 : public Analysis {

  public:

    /// Default constructor.
    D0_2007_S7075677();


    /// Factory method 
    static Analysis* create() {
      return new D0_2007_S7075677();
    }


    /// @name Analysis methods
    //@{ 
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D * _h_yZ;
    //@}

  };


}

#endif
