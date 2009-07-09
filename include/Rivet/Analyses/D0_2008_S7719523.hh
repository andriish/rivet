// -*- C++ -*-
#ifndef RIVET_D0_2008_S7719523_HH
#define RIVET_D0_2008_S7719523_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief Measurement of isolated gamma + jet + X differential cross-sections
  /// Inclusive isolated gamma + jet cross-sections, differential in pT(gamma), for 
  /// various photon and jet rapidity bins.
  ///
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  class D0_2008_S7719523 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
     D0_2008_S7719523();

    /// Factory method 
    static Analysis* create() {
      return new D0_2008_S7719523();
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
    AIDA::IHistogram1D* _h_central_same_cross_section;
    AIDA::IHistogram1D* _h_central_opp_cross_section;
    AIDA::IHistogram1D* _h_forward_same_cross_section;
    AIDA::IHistogram1D* _h_forward_opp_cross_section;
    //@}

  };

}

#endif
