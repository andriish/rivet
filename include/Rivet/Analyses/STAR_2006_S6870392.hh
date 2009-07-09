// -*- C++ -*-
#ifndef RIVET_STAR_2006_S6870392_HH
#define RIVET_STAR_2006_S6870392_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief inclusive jet cross-section in pp at 200 GeV
  class STAR_2006_S6870392 : public Analysis {

  public:

    /// @name Construction
    //@{
    /// Constructor
    STAR_2006_S6870392();

    /// Factory method 
    static Analysis* create() {
      return new STAR_2006_S6870392();
    }
    //@}


    /// @name Analysis methods
    //@{ 
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}


  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D * _h_jet_pT_MB;
    AIDA::IHistogram1D * _h_jet_pT_HT;
    //@}

  };


}

#endif
