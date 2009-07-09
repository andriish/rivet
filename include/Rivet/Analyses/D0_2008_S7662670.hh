// -*- C++ -*-
#ifndef RIVET_D0_2008_S7662670_HH
#define RIVET_D0_2008_S7662670_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief Measurement of D0 differential jet cross sections
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  class D0_2008_S7662670 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    D0_2008_S7662670();

    /// Factory method 
    static Analysis* create() {
      return new D0_2008_S7662670();
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
    AIDA::IHistogram1D* _h_dsigdptdy_y00_04;
    AIDA::IHistogram1D* _h_dsigdptdy_y04_08;
    AIDA::IHistogram1D* _h_dsigdptdy_y08_12;
    AIDA::IHistogram1D* _h_dsigdptdy_y12_16;
    AIDA::IHistogram1D* _h_dsigdptdy_y16_20;
    AIDA::IHistogram1D* _h_dsigdptdy_y20_24;
    //@}

  };

}

#endif
