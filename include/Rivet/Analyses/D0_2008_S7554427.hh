// -*- C++ -*-
#ifndef RIVET_D0_2008_S7554427_HH
#define RIVET_D0_2008_S7554427_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief Measurement of D0 Run II Z pT differential cross-section shape
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  /// @author Frank Siegert
  class D0_2008_S7554427 : public Analysis {

  public:

    /// Default constructor.
    D0_2008_S7554427();


    /// Factory method 
    static Analysis* create() {
      return new D0_2008_S7554427();
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
    AIDA::IHistogram1D * _h_ZpT;
    AIDA::IHistogram1D * _h_forward_ZpT;
    //@}

  };


}

#endif
