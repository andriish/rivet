// -*- C++ -*-
#ifndef RIVET_STAR_2008_S7993412_HH
#define RIVET_STAR_2008_S7993412_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief di-hadron correlations in d-Au at 200 GeV
  class STAR_2008_S7993412 : public Analysis {

  public:

    /// @name Construction
    //@{
    /// Constructor
    STAR_2008_S7993412();

    /// Factory method 
    static Analysis* create() {
      return new STAR_2008_S7993412();
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
    AIDA::IProfile1D * _h_Y_jet_trigger;
    AIDA::IProfile1D * _h_Y_jet_associated;
    //@}

  };


}

#endif
