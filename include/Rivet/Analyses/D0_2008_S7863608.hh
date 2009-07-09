// -*- C++ -*-
#ifndef RIVET_D0_2008_S7863608_HH
#define RIVET_D0_2008_S7863608_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief Measurement differential Z/gamma* + jet +X cross sections
  /// @author Gavin Hesketh, Andy Buckley, Frank Siegert
  class D0_2008_S7863608 : public Analysis {

  public:

    /// @name Construction
    //@{
    /// Constructor
    D0_2008_S7863608();

    /// Factory method 
    static Analysis* create() {
      return new D0_2008_S7863608();
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
    AIDA::IHistogram1D * _h_jet_pT_cross_section;
    AIDA::IHistogram1D * _h_jet_y_cross_section;
    AIDA::IHistogram1D * _h_Z_pT_cross_section;
    AIDA::IHistogram1D * _h_Z_y_cross_section;
    AIDA::IHistogram1D * _h_total_cross_section;
    //@}

  };


}

#endif
